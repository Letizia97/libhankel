.. _tabulated-form-factors:

Tabulated form factors (for C users)
======================================

When the form factor is a set of points rather than a formula -- measured data,
or the output of another code -- :c:func:`tabulated_ff_create` turns it into
something ``hankel_transform`` can call. Hand it the table and a rule for the
high-:math:`q` tail and it returns a handle you pass straight through as
``f_ctx``; see :ref:`tabulated_ff <tabulated_ff_c_api>` for the API and
:ref:`c-examples-tabulated-f` for a complete program.

Interpolating is the easy part. The part that needs care is what the form
factor does *outside* the tabulated range, and this page is mostly about that.


Why the tail matters
----------------------

``hankel_transform`` does not evaluate the form factor at the points you pass
in ``x``. It evaluates it at quadrature nodes chosen by the strategy, and those
go a very long way outside any real dataset:

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

No experimental dataset covers 35 decades, so the callback is asked for values
outside the table on every single call. :c:func:`cubic_interp_eval` returns NaN
there, because Boost refuses to extrapolate.

.. warning::

    Feeding a bare spline to ``hankel_transform`` fails **silently**. The
    status code is 0 in every row below.

Transforming a table of the ``g_dab`` form factor sampled over
:math:`q \in [10^{-3}, 10]`, with ``n_eval = 250``, ``eps_rel = 1e-9`` and
``f_max = 0.1``:

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

The NaN rows are obvious enough once you look at the output. The last two are
the dangerous ones: the extrapolation-driven NaNs get absorbed by the
convergence machinery, so the call returns finite, plausible-looking numbers
that are wrong by about ten percent.


Choosing the tail
-------------------

Above the table you choose one of ``tabulated_tail_t``. The integrand is
:math:`f(q)\,J_0(qr)\,q`, and for large argument

.. math::

    J_0(z) \sim \sqrt{\frac{2}{\pi z}} \cos\left(z - \frac{\pi}{4}\right),

so a tail :math:`f(q) \sim q^{-p}` leaves an integrand going like
:math:`q^{1/2 - p}`. That gives three possibilities, of which the library
offers the two that work:

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

:c:func:`tabulated_ff_create` enforces the threshold: an exponent at or below
3/2 is rejected with -14 rather than turned into a wrong answer later. Note
that clamping to the endpoint value, the intuitive thing to reach for, is
:math:`p = 0` -- the one option that is outright broken.

Porod's law gives :math:`p = 4` for a sharp interface. If you would rather not
assume that, pass ``0.0`` as the exponent and it is fitted from the last points
of the table by least squares in log-log space.

.. note::

    A fit can only see the slope it was shown. The example's table stops at
    :math:`q = 1` and fits :math:`p = 3.82`; the same curve carried out to
    :math:`q = 10` fits :math:`3.997`, against an analytic value of exactly 4.
    So a -14 from a fitted exponent usually means the table stops before the
    asymptotic regime, not that the physics diverges -- extend the table, or
    pass the exponent you know applies.

Below the table the value is held flat at the first tabulated point, the
:math:`q \to 0` plateau. This is not configurable and does not need to be: the
integrand's factor of :math:`q` suppresses that region regardless.


Cost
------

Each evaluation does a binary search plus a cubic evaluation instead of a
handful of flops:

.. table:: Cost per form factor evaluation.

    +------------------------------------------+-------------------+
    | Form factor                              | Per evaluation    |
    +==========================================+===================+
    | Analytic built-in (sphere)               | ~10 ns            |
    +------------------------------------------+-------------------+
    | PCHIP spline, 1000-point table           | ~260 ns           |
    +------------------------------------------+-------------------+

The cost grows slowly with table size (roughly 55 ns per decade of points,
which is the binary search), so a coarser table helps only marginally. Building
the handle is a one-off ~38 microseconds for a 1000-point table: build it once
outside the loop and reuse it, never inside the callback. For a strategy making
tens of thousands of evaluations this puts the transform in the
several-millisecond range rather than the sub-millisecond range.
