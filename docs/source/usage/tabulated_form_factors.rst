.. _tabulated-form-factors:

Tabulated form factors (for C users)
======================================

When the form factor is a set of discrete points rather than a formula — measured
data, or the output of another simulation — :c:func:`tabulated_ff_create` wraps it
into a callable that ``hankel_transform`` can use. You supply two things: the table
itself, and a rule for the high-:math:`q` tail beyond the last tabulated point.
The function returns a handle you pass as ``f_ctx``; see
:ref:`tabulated_ff <tabulated_ff_c_api>` for the API and
:ref:`c-examples-tabulated-f` for a complete worked example.

Using the tabulated form factor involves two decisions: which **interpolation
method** to use within the table, and which **tail rule** to apply outside it.
Interpolation is straightforward; the tail requires more care, and the rest of
this page is largely about that.


Choosing an interpolation method
---------------------------------

Within the tabulated range, the library offers three interpolators to choose from
when building the form factor with :c:func:`tabulated_ff_create`. Each makes
different trade-offs:

.. list-table:: Interpolation options.
   :header-rows: 1

   * - Method
     - Implementation
     - Error
     - Best for
   * - Linear
     - :c:func:`linear_interp_create`
     - O(h²)
     - Speed-critical code or small tables
   * - Cubic (PCHIP)
     - :c:func:`cubic_interp_create`
     - O(h⁴)
     - Smooth curves or when kinks matter
   * - Log-linear
     - :c:func:`loglinear_interp_create`
     - Constant relative
     - Form factors across wide ranges (SAX/SANS)

**Linear** is simple and fast. The interpolated values are guaranteed to stay
within the bracket of the tabulated endpoints, and monotonicity is preserved:
non-negative data in, non-negative data out. The price is a discontinuous slope
at every tabulated point.

**Cubic (PCHIP)** produces smooth curves through the data with no overshoot or
undershoot, at a small cost in evaluation time. Use this when the shape of the
interpolated curve between tabulations matters, or when feeding the form factor
to visualization or plotting code.

**Log-linear** interpolates in :math:`\log(y)` space, returning :math:`\exp(\text{result})`.
This gives exponential (power-law) interpolation, which is ideal for form factors
spanning orders of magnitude: it provides constant *relative* error across
decades rather than constant absolute error. For a form factor dropping from
1 to :math:`10^{-10}`, linear gives 1e-10 absolute error everywhere (tiny at
low q, wasted at high q); log-linear gives the same relative error at every
scale. **This is the recommended choice for SAXS and SANS form factors.**
Log-linear requires all tabulated values to be strictly positive.

See the :ref:`Interpolation section of the C API <interpolation_c_api>` for detailed documentation on each interpolator.


Choosing the tail
-------------------

Beyond the last tabulated point, you select a tail rule via ``tabulated_tail_t``.
The integrand is :math:`f(q)\,J_0(qr)\,q`, and for large argument

.. math::

    J_0(z) \sim \sqrt{\frac{2}{\pi z}} \cos\left(z - \frac{\pi}{4}\right),

so a tail :math:`f(q) \sim q^{-p}` leaves an integrand behaving like
:math:`q^{1/2 - p}`. This gives three convergence cases, of which the library
offers the two that are well-behaved:

.. table:: High-q tail options.

    +------------------------------+---------+----------------------------------------------+
    | Tail                         | p       | Behaviour                                    |
    +==============================+=========+==============================================+
    | Hold flat at f(x[n-1])       | 0       | **Does not converge.** Not offered: the      |
    |                              |         | integrand oscillates with amplitude growing  |
    |                              |         | like sqrt(q).                                |
    +------------------------------+---------+----------------------------------------------+
    | ``TABULATED_TAIL_ZERO``      | --      | Converges, but the discontinuity at the      |
    |                              |         | cutoff puts ringing into the result. Right   |
    |                              |         | only when the data really is zero past the   |
    |                              |         | table, e.g. a hard instrumental cutoff.      |
    +------------------------------+---------+----------------------------------------------+
    | ``TABULATED_TAIL_POWER_LAW`` | > 3/2   | Converges, and continuous at the join.       |
    |                              |         | Recommended.                                 |
    +------------------------------+---------+----------------------------------------------+

:c:func:`tabulated_ff_create` enforces this threshold: an exponent at or below
3/2 is rejected with an error code (-14) rather than silently producing a wrong
answer. Note that clamping to the endpoint value — the intuitive fallback —
corresponds to :math:`p = 0`, which is the one option that does not converge.

Porod's law gives :math:`p = 4` for a sharp interface. If that assumption does
not apply, pass ``0.0`` as the exponent and the library will fit it from the
last few points of the table by least squares in log-log space.

.. note::

    A fit can only see the slope it was shown. The example's table stops at
    :math:`q = 1` and fits :math:`p = 3.82`; the same curve extended to
    :math:`q = 10` fits :math:`3.997`, against an analytic value of exactly 4.
    An error code (-14) from a fitted exponent therefore usually means the table
    stops before the asymptotic regime, not that the physics diverges — either
    extend the table further into the tail, or pass the exponent you know applies.

Below the table the value is held flat at the first tabulated point, the
:math:`q \to 0` plateau. This is not configurable and does not need to be: the
integrand's factor of :math:`q` suppresses that region regardless.


Why ``tabulated_ff`` is needed
------------------------------

``hankel_transform`` does not evaluate the form factor at the points you pass in
``x``. It evaluates it at quadrature nodes chosen by the strategy, and those reach
far outside the range of any real dataset:

.. table:: Arguments requested by each strategy, over 19 output points with ``n_eval = 250``.

    +----------------------+----------------------------+
    | Strategy             | Range of q requested       |
    +======================+============================+
    | DHT_Anderson_801     | 1.1e-15  to  3.3e+20       |
    +----------------------+----------------------------+
    | DHT_Key_101          | 7.0e-06  to  1.2e+02       |
    +----------------------+----------------------------+
    | DHT_Guptasarma       | 5.2e-11  to  1.6e+01       |
    +----------------------+----------------------------+
    | Fixed_DE_Ogata       | 3.7e-04  to  5.2e+01       |
    +----------------------+----------------------------+
    | Adaptive_DE_Ooura    | 1.4e-08  to  7.1e+00       |
    +----------------------+----------------------------+
    | QWE_Key              | 9.1e-11  to  2.3e+00       |
    +----------------------+----------------------------+
    | QWE_Chave            | 3.7e-10  to  2.3e+00       |
    +----------------------+----------------------------+

No experimental dataset covers 35 decades. This means a bare interpolator —
one with no tail rule — will return NaN for out-of-range queries on essentially
every call. The consequences are not always obvious: some strategies propagate
NaN to the output immediately, but others absorb it silently inside their
convergence loop and return a finite, plausible-looking number that is simply
wrong.

The table below shows this concretely, using a table of the ``g_dab`` form
factor sampled over :math:`q \in [10^{-3}, 10]`, with ``n_eval = 250``,
``eps_rel = 1e-9`` and ``f_max = 0.1``. The status code is 0 in every case —
the bare spline never raises an error.

.. table:: Bare spline vs. the same table through ``tabulated_ff``.

    +----------------------+-------------------------+---------------------------+
    | Strategy             | Bare spline             | Through ``tabulated_ff``  |
    +======================+=========================+===========================+
    | DHT_Key_101          | all 19 points NaN       | max rel. error 3.4e-04    |
    +----------------------+-------------------------+---------------------------+
    | Fixed_DE_Ogata       | 18 of 19 points NaN     | max rel. error 4.0e-04    |
    +----------------------+-------------------------+---------------------------+
    | Adaptive_DE_Ooura    | all 19 points NaN       | max rel. error 1.4e-04    |
    +----------------------+-------------------------+---------------------------+
    | QWE_Key              | no NaN, **9.8 % wrong** | max rel. error 1.9e-05    |
    +----------------------+-------------------------+---------------------------+
    | QWE_Chave            | no NaN, **9.8 % wrong** | max rel. error 1.9e-05    |
    +----------------------+-------------------------+---------------------------+

The NaN rows are at least detectable once you inspect the output. The QWE rows
are the dangerous case: the NaNs are absorbed by the convergence machinery and
the call returns finite numbers that are wrong by ten percent with no indication
that anything went wrong.
