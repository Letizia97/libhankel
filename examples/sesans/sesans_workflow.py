"""
The SESANS workflow, and where the Hankel transform sits inside it.

A SESANS (spin-echo small-angle neutron scattering) measurement, e.g. on Larmor
at the ISIS Neutron and Muon Source, does not measure I(q) directly: it measures
the polarisation of the neutron beam as a function of the spin-echo length xi.
The model side of the analysis therefore has to turn a scattering model I(q)
into a real-space correlation function G(xi), and that step is a zeroth-order
Hankel transform:

    G(xi) = 1 / (2 pi) * integral_0^inf q J0(q xi) I(q) dq
    P(xi) = G(xi) - G(0)

In SasView this is done by ``sasmodels.sesans.SesansTransform``, which builds a
dense matrix H of the Hankel kernel on a log-spaced q grid and applies it with a
dot product. The call chain, starting from anything that fits or plots SESANS
data, is:

    sasmodels.direct_model.DirectModel(data, model)     # data.isSesans is True
      -> DataMixin._interpret_data
        -> direct_model._make_sesans_transform(data)
          -> sesans.SesansTransform(z, SElength, lam, zaccept, Rmax)
            -> SesansTransform._set_hankel   # builds H and H0
    calculator(**pars)
      -> DataMixin._calc_theory
        -> call_kernel(...)                  # I(q) on transform.q_calc
        -> SesansTransform.apply(Iq)         # G = H.T @ Iq, G0 = H0 @ Iq

This script runs that workflow for a dilute dispersion of spheres, then computes
the same P(xi) with LibHankel doing the transform instead, and compares the two.

Requires ``sasmodels`` alongside LibHankel:

    python -m pip install sasmodels
    python examples/sesans/sesans_workflow.py
"""

import libhankel
import numpy as np
from sasmodels.core import load_model
from sasmodels.data import empty_sesans
from sasmodels.direct_model import DirectModel

# Sample: a dilute dispersion of hard spheres.
RADIUS = 1000.0  # A
CONTRAST = 1.0  # sasmodels: 1e-6 A^-2; libhankel: arbitrary units (see below)
WAVELENGTH = 5.0  # A

# Spin-echo lengths, in A. A sphere correlation function dies at xi = 2 R,
# so this range covers the whole of G(xi) for the radius above.
SPIN_ECHO_LENGTHS = np.logspace(np.log10(50.0), np.log10(5000.0), 50)

# LibHankel strategy used in place of the SesansTransform matrix.
STRATEGY = "QWE_Chave"
STRATEGY_PARAMS = {"n_eval": 250, "eps_rel": 1e-9}


def sesans_with_sasmodels(xi):
    """
    Run the sasmodels SESANS workflow, the one that calls SesansTransform.

    Returns P(xi) = G(xi) - G(0), G(0), and the q grid the transform used.
    """

    # SesansData carries isSesans = True, which is what makes DirectModel
    # build a SesansTransform rather than a resolution function.
    data = empty_sesans(z=xi, wavelength=WAVELENGTH)
    model = load_model("sphere")
    calculator = DirectModel(data, model)

    polarisation = calculator(
        radius=RADIUS,
        sld=CONTRAST,
        sld_solvent=0.0,
        scale=1.0,
        background=0.0,
    )

    # calculator.resolution is the SesansTransform built above; _H0 is the
    # xi = 0 row of the Hankel kernel, so this is the G(0) it subtracted.
    transform = calculator.resolution
    g_zero = np.dot(transform._H0, calculator.Iq_calc)

    return polarisation, g_zero, transform.q_calc


def sesans_with_libhankel(xi):
    """
    Same calculation, with LibHankel evaluating the Hankel transform.

    No q grid and no matrix: LibHankel integrates q I(q) J0(q xi) over the whole
    half-line for each xi, using the built-in "sphere" form factor

        I(q) = [4 pi eta (sin(qR) - qR cos(qR)) / q^3]^2

    which is the same I(q) as the sasmodels sphere up to a constant factor.
    """

    params = np.array([RADIUS, CONTRAST])

    g_xi = np.asarray(
        libhankel.hankel_transform(0, "sphere", xi, params, STRATEGY, STRATEGY_PARAMS)
    ) / (2 * np.pi)

    # J0(0) = 1, so G(0) is just the integral of q I(q), which for a sphere is
    # 2 pi eta^2 R^4 exactly. For a form factor with no closed form, integrate
    # q I(q) numerically instead (e.g. scipy.integrate.quad, split over the
    # oscillations of the integrand).
    g_zero = 2 * np.pi * CONTRAST**2 * RADIUS**4

    return g_xi - g_zero, g_zero


def main():
    xi = SPIN_ECHO_LENGTHS

    p_sasmodels, g0_sasmodels, q_calc = sesans_with_sasmodels(xi)
    p_libhankel, g0_libhankel = sesans_with_libhankel(xi)

    print(f"sphere of radius {RADIUS:g} A, {len(xi)} spin-echo lengths")
    print(
        f"sasmodels : SesansTransform on {len(q_calc):d} q points, "
        f"q = {q_calc[0]:.3g} to {q_calc[-1]:.3g} 1/A"
    )
    print(f"libhankel : {STRATEGY} over the full q range, no grid")

    # The two I(q) differ by a constant factor (sasmodels normalises by volume
    # and converts to cm^-1), so compare the normalised correlation function
    # P(xi) / G(0), which is independent of that factor.
    norm_sasmodels = p_sasmodels / g0_sasmodels
    norm_libhankel = p_libhankel / g0_libhankel
    difference = norm_libhankel - norm_sasmodels

    print()
    print(f"{'xi (A)':>10} {'sasmodels':>14} {'libhankel':>14} {'difference':>12}")
    for value, sas, lib, diff in zip(xi, norm_sasmodels, norm_libhankel, difference):
        print(f"{value:10.1f} {sas:14.6f} {lib:14.6f} {diff:12.2e}")

    print()
    print(f"largest absolute difference in P(xi)/G(0): {np.max(np.abs(difference)):.2e}")


if __name__ == "__main__":
    main()
