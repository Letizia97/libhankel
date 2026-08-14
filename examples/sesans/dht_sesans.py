"""
A drop-in replacement for sasmodels' ``SesansTransform``, built on a LibHankel
digital filter.

``sasmodels.sesans.SesansTransform`` computes

.. math::   G(\\xi) = \\frac{1}{2\\pi} \\int_0^\\infty q\\, J_0(q\\xi)\\, I(q)\\, dq

by brute force: it lays down a dense log-spaced grid of q (43936 points for the
ISIS dataset), precomputes a matrix of :math:`J_0` values on it, and reduces
each model evaluation to a matrix product. The matrix is 20 MB and takes about
50 ms to build.

A digital filter computes the same integral as a weighted sum over abscissae
that are a fixed table divided by :math:`\\xi`::

    lambda = KK201Hankel[i][0] / x;
    res = res + (*f)(lambda, f_ctx) * lambda * KK201Hankel[i][ind] / x;

(``src/hankel_DHT.c``). The abscissae depend only on the filter and on
:math:`\\xi`, so they are known before the model is evaluated -- which is
exactly the contract ``q_calc`` expresses. That makes a filter, alone among the
LibHankel strategies, expressible as a sasmodels resolution object: the
adaptive strategies choose their abscissae from values they have already seen,
so they cannot say in advance where the model is needed.

This module reads the abscissae and weights through :func:`libhankel.dht_nodes`
and does the weighted sum in numpy. Nothing crosses into C per model
evaluation, so there is no callback cost; LibHankel supplies the quadrature
rule and numpy applies it.

Run ``benchmark_dht_transform.py`` for timings and accuracy against the real
``.ses`` files.
"""

import libhankel
import numpy as np
from numpy import pi

#: Filters whose abscissae are tabulated, longest last. Longer filters cost
#: proportionally more model evaluations.
FILTERS = ("DHT_Key_51", "DHT_Key_101", "DHT_Key_201", "DHT_Anderson_801")


class DHTSesansTransform:
    """
    Spin-Echo SANS transform calculator, using a LibHankel digital filter.

    Interchangeable with :class:`sasmodels.sesans.SesansTransform`: it exposes
    the same ``q``, ``q_calc`` and ``apply`` that
    ``DirectModel._calc_theory`` uses, and takes the same five positional
    arguments that ``direct_model._make_sesans_transform`` passes.

    *z* is the spin-echo length in the original data units. The transform does
    not use it; the SasView GUI does, so it is stored unchanged.

    *SElength* (A) is the set of spin-echo lengths in the measured data.

    *lam* (A) is the wavelength, either a scalar or one value per spin-echo
    length (time-of-flight instruments vary it point by point).

    *zaccept* (1/A) is the maximum acceptance of the scattering vector in the
    spin-echo encoding direction.

    *Rmax* (A) is accepted for signature compatibility and ignored, as it is by
    ``SesansTransform._set_hankel``. A filter has no grid to size.

    *strategy* names the filter, one of :data:`FILTERS`.
    """

    #: SElength in the original data units; not used by the transform, but the
    #: GUI reads it, so make sure that it is present.
    q = None  # type: np.ndarray

    #: q values to calculate when computing the transform.
    q_calc = None  # type: np.ndarray

    def __init__(self, z, SElength, lam, zaccept, Rmax, strategy="DHT_Key_201"):
        # type: (np.ndarray, np.ndarray, np.ndarray, float, float, str) -> None
        self.q = z
        self.strategy = strategy
        self._set_filter(SElength, lam, zaccept)

    def apply(self, Iq):
        # type: (np.ndarray) -> np.ndarray
        """
        Apply the SESANS transform to the computed I(q).

        ``SesansTransform.apply`` writes both reductions as ``np.dot``. That is
        right for its (n_q, n_xi) matrix product, but for a vector reduction it
        is a trap: numpy sends ``np.dot`` to BLAS, and on an array this size
        threaded BLAS spends far longer starting and joining threads than it
        does on the arithmetic. Measured here, ``np.dot(H0, Iq)`` on 11457
        points costs 2.9 ms against 0.015 ms for the multiply and sum below --
        and ``np.dot`` of two vectors of ones is just as slow, so it is
        overhead, not the values.
        """
        G0 = np.sum(self._H0 * Iq)
        G = np.sum(np.reshape(Iq, self._H.shape) * self._H, axis=0)
        P = G - G0
        return P

    def _set_filter(self, SElength, lam, zaccept):
        # type: (np.ndarray, np.ndarray, float) -> None
        SElength = np.asarray(SElength, dtype=float)
        lam = np.broadcast_to(np.asarray(lam, dtype=float), SElength.shape)

        abscissae, weights = libhankel.dht_nodes(self.strategy, 0)
        abscissae = np.asarray(abscissae)
        weights = np.asarray(weights)

        # The filter needs I(q) at node/xi, one column of q per spin-echo
        # length. Both the abscissae and the weights carry a 1/xi.
        q = abscissae[:, None] / SElength[None, :]
        H = q * weights[:, None] / SElength[None, :] / (2 * pi)

        # Same acceptance masking as SesansTransform._set_hankel: neutrons
        # scattered beyond the analyser acceptance never reach the detector, so
        # they contribute nothing. When q lam / 2 pi > 1 the wavelength is too
        # large to reach that q at any angle, which gives theta = NaN; reverse
        # the condition rather than testing theta > zaccept, because every
        # comparison with NaN is False.
        with np.errstate(invalid="ignore"):
            theta = np.arcsin(q * lam[None, :] / (2 * pi))
        H[~(theta <= zaccept)] = 0

        self.q_calc = q.ravel()
        self._H = H
        self._H0 = self._g0_weights(self.q_calc)

    @staticmethod
    def _g0_weights(q_calc):
        # type: (np.ndarray) -> np.ndarray
        """
        Quadrature weights for G(0) over the abscissae the filter already uses.

        ``apply`` returns ``G(xi) - G(0)``, and G(0) is the limit of the
        transform as xi goes to zero. A filter cannot be evaluated there --
        its abscissae are ``node/xi`` -- so G(0) needs its own quadrature.
        ``SesansTransform`` gets it free from its dense grid
        (``H0 = dq/(2 pi) * q``, a rectangle rule); a filter has no such grid.

        But G(0) is a plain integral, not a Hankel transform:

        .. math::   G(0) = \\frac{1}{2\\pi} \\int_0^\\infty q\\, I(q)\\, dq
                         = \\frac{1}{2\\pi} \\int q^2 I(q)\\, d\\ln q

        and the union of the abscissae over all the spin-echo lengths already
        covers the q range densely -- 11457 points across 13 decades for the
        ISIS dataset, because each spin-echo length shifts the same node table
        by a different 1/xi. So a trapezoid rule in log q over that union costs
        no extra model evaluations at all, and on the ISIS data it is more
        accurate than the rectangle rule it replaces.

        The one case this does not cover is a single spin-echo length, where
        the union is just the node table itself and the sampling is roughly 19
        points per decade.

        Like ``H0`` in ``SesansTransform``, these weights are not acceptance
        masked.
        """
        order = np.argsort(q_calc)
        log_q = np.log(q_calc[order])

        # Trapezoid weights: each point takes half of the interval on either
        # side. Repeated abscissae give a zero-width interval and drop out.
        step = np.diff(log_q)
        trapezoid = np.zeros_like(log_q)
        trapezoid[:-1] += step / 2
        trapezoid[1:] += step / 2

        weights = np.empty_like(q_calc)
        weights[order] = trapezoid * q_calc[order] ** 2 / (2 * pi)
        return weights


def use_dht_transform(strategy="DHT_Key_201"):
    """
    Patch sasmodels so that every SESANS ``DirectModel`` uses a filter.

    ``DirectModel`` builds its transform in
    ``direct_model._make_sesans_transform``, which reads the spin-echo lengths,
    wavelengths and analyser acceptance off the data and converts them to
    Angstroms and radians. Only the final line picks the transform, so
    replacing the function keeps all of that unit handling and swaps the
    quadrature underneath it.

    Returns the original function, so a caller can restore it.
    """
    from sasmodels import direct_model

    original = direct_model._make_sesans_transform

    def _make_dht_transform(data):
        SElength, SEunits = data.x, data._xunit
        wavelength, wunits = data.source.wavelength, data.source.wavelength_unit
        theta_max, theta_units = data.sample.zacceptance
        if SEunits != "A" or wunits != "A" or theta_units != "radians":
            from sasdata.data_util.nxsunit import Converter

            SElength = Converter("A")(SElength, units=SEunits)
            wavelength = Converter("A")(wavelength, units=wunits)
            theta_max = Converter("radian")(theta_max, units=theta_units)

        Rmax = 10000000
        zaccept = 2 * np.pi / np.max(wavelength) * np.sin(theta_max)
        return DHTSesansTransform(data.x, SElength, wavelength, zaccept, Rmax, strategy=strategy)

    direct_model._make_sesans_transform = _make_dht_transform
    return original
