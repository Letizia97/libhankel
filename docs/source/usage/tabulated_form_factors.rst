.. _tabulated-form-factors:

Tabulated form factors (for C users)
======================================

Sometimes the form factor is not a formula but a set of points -- measured
data, or the output of another code. LibHankel provides a PCHIP cubic
interpolator (see :ref:`interpolation_c_api`) to turn such a table into
something ``hankel_transform`` can call.

Interpolating is the easy part. The part that needs care is deciding what the
form factor does *outside* the tabulated range, and this page is mostly about
that.


Why the tails matter
----------------------

``hankel_transform`` does not evaluate the form factor at the points you pass
in ``x``. It evaluates it at quadrature nodes chosen by the strategy, and
those nodes go a very long way outside any real dataset. Measured over the
example in :ref:`c-examples-tabulated-f` (19 output points, ``n_eval = 250``),
the strategies asked for these ranges of :math:`q`:

.. table:: Range of form factor arguments requested by each strategy.

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

No experimental dataset covers 35 decades. So the form factor callback will be
asked for values outside the table, on every single call, no matter how wide
the table is.


What goes wrong if you ignore them
------------------------------------

:c:func:`cubic_interp_eval` returns NaN outside ``[x[0], x[n-1]]``, because
Boost refuses to extrapolate. Passing the bare interpolator to
``hankel_transform`` therefore gives one of two bad outcomes, both of which
still return status code 0:

.. warning::

    Feeding an un-extended spline to ``hankel_transform`` fails **silently**.
    The status code is 0 in every case below.

Transforming a table of the ``g_dab`` form factor sampled over
:math:`q \in [10^{-3}, 10]`, first with the bare spline and then with the
tails defined as described further down:

.. table:: Bare spline vs. spline with defined tails.

    +----------------------+-------------------------+---------------------------+
    | Strategy             | Bare spline             | With tails defined        |
    +======================+=========================+===========================+
    | DHT_Key_101          | all 19 points NaN       | max rel. error 3.4e-04    |
    +----------------------+-------------------------+---------------------------+
    | Fixed_DE_Ogata       | all 19 points NaN       | max rel. error 6.2e-04    |
    +----------------------+-------------------------+---------------------------+
    | Adaptive_DE_Ooura    | all 19 points NaN       | max rel. error 1.4e-04    |
    +----------------------+-------------------------+---------------------------+
    | QWE_Chave            | no NaN, **9.8 % wrong** | max rel. error 1.9e-05    |
    +----------------------+-------------------------+---------------------------+

The NaN rows are obvious enough once you look at the output. The ``QWE_Chave``
row is the dangerous one: the extrapolation-driven NaNs get absorbed by the
convergence machinery, so the call returns finite, plausible-looking numbers
that are wrong by about ten percent.


Choosing the tails
--------------------

Below the table
^^^^^^^^^^^^^^^^^

Hold the value flat at ``f(x[0])``. This is the :math:`q \to 0` plateau, and
it is numerically harmless: the Hankel integrand carries a factor of :math:`q`
which suppresses that region regardless.

Above the table
^^^^^^^^^^^^^^^^^

This one is a real choice, and it affects whether the integral converges at
all. The integrand is :math:`f(q)\,J_0(qr)\,q`, and for large argument

.. math::

    J_0(z) \sim \sqrt{\frac{2}{\pi z}} \cos\left(z - \frac{\pi}{4}\right),

so an assumed high-:math:`q` behaviour of :math:`f(q) \sim q^{-p}` leaves an
integrand going like :math:`q^{1/2 - p}`. That gives three options:

.. table:: High-q tail options.

    +--------------------------+---------+--------------------------------------------------+
    | Tail                     | p       | Behaviour                                        |
    +==========================+=========+==================================================+
    | Hold flat at f(x[n-1])   | 0       | **Does not converge.** The integrand oscillates  |
    |                          |         | with amplitude growing like sqrt(q).             |
    +--------------------------+---------+--------------------------------------------------+
    | Truncate to zero         | --      | Converges, but the discontinuity at the cutoff   |
    |                          |         | puts ringing into the result.                    |
    +--------------------------+---------+--------------------------------------------------+
    | Power law f ~ q^-p       | > 3/2   | Converges, and continuous at the join.           |
    |                          |         | Recommended.                                     |
    +--------------------------+---------+--------------------------------------------------+

Porod's law gives :math:`p = 4` for a sharp interface, comfortably inside the
convergent range. If you do not want to assume that, fit the slope from the
last part of the table in log-log space, which is what
:ref:`c-examples-tabulated-f` does -- it recovers :math:`p = 3.996` against an
analytic value of exactly 4.

Note that clamping to the endpoint value, which is the intuitive thing to
reach for, is the one option that is outright broken.


Cost
------

An interpolated form factor is slower per evaluation than an analytic one,
because each call does a binary search plus a cubic evaluation instead of a
handful of flops. Measured on this library:

.. table:: Cost per form factor evaluation.

    +------------------------------------------+-------------------+
    | Form factor                              | Per evaluation    |
    +==========================================+===================+
    | Analytic built-in (sphere)               | ~10 ns            |
    +------------------------------------------+-------------------+
    | PCHIP spline, 1000-point table           | ~260 ns           |
    +------------------------------------------+-------------------+

The spline evaluation cost grows slowly with table size (roughly 55 ns per
decade of points, which is the binary search). Building the spline is a
one-off cost of about 38 microseconds for a 1000-point table, so build it once
outside the loop and reuse the handle -- never rebuild it inside the callback.

For a strategy that makes tens of thousands of evaluations this puts the
transform in the several-millisecond range rather than the sub-millisecond
range. That is usually acceptable; if it is not, a coarser table helps only
marginally, since most of the cost is a fixed per-call overhead rather than
the search.


Worked example
----------------

See :ref:`c-examples-tabulated-f` for a complete, self-checking program that
builds a spline from a table, defines both tails, transforms it, and compares
the answer against the analytic transform of the same form factor.
