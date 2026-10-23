.. _strategy-selection:

Strategy Selection
================================

LibHankel provides several implementations of the Hankel transform. All of them are
accessed through the same ``hankel_transform`` function by passing a ``strategy_name``
string, so switching strategies is a one-word change with no other code to modify.

**Not sure which to pick?** Start with ``FFTLog``. It requires no strategy parameters,
runs in O(N log N) time across the whole output array, and is typically the most accurate
option for smooth, slowly-varying form factors evaluated over a wide dynamic range.
If your form factor is strongly oscillatory or has sharp features, try one of the adaptive
strategies (``QWE_Chave``, ``QWE_Key``, ``Adaptive_DE_Ooura``) instead — they let you
control convergence via a tolerance parameter.

The following table lists all available strategies.
*Nodes* records whether the quadrature points are fixed or chosen adaptively (see
:ref:`fixed-vs-adaptive` below).
The *Index* column is the numbering used in SASfit, for cross-reference.


.. table:: Hankel transform strategies available in LibHankel.

    +----------------------+---------------------------------------+----------+----------------------+
    | Strategy Name        | Type / Method                         | Nodes    | Index (as in SASfit) |
    +======================+=======================================+==========+======================+
    | DHT_Guptasarma       | Digital filter (Guptasarma & Singh)   | fixed    | 6                    |
    +----------------------+---------------------------------------+----------+----------------------+
    | DHT_Guptasarma_Fast  | Digital filter (Guptasarma & Singh)   | fixed    | 7                    |
    +----------------------+---------------------------------------+----------+----------------------+
    | DHT_Key_51           | Digital filter (Key)                  | fixed    | 8                    |
    +----------------------+---------------------------------------+----------+----------------------+
    | DHT_Key_101          | Digital filter (Key)                  | fixed    | 9                    |
    +----------------------+---------------------------------------+----------+----------------------+
    | DHT_Key_201          | Digital filter (Key)                  | fixed    | 10                   |
    +----------------------+---------------------------------------+----------+----------------------+
    | DHT_Anderson_801     | Digital filter (Anderson)             | fixed    | 11                   |
    +----------------------+---------------------------------------+----------+----------------------+
    | Fixed_DE_Ogata       | Double-exponential quadrature         | fixed    | 1                    |
    |                      | at the zeros of the Bessel function   |          |                      |
    +----------------------+---------------------------------------+----------+----------------------+
    | Adaptive_DE_Ooura    | Double-exponential quadrature         | adaptive | 0                    |
    +----------------------+---------------------------------------+----------+----------------------+
    | QWE_Key              | QWE with continued fraction expansion | adaptive | 12                   |
    +----------------------+---------------------------------------+----------+----------------------+
    | QWE_Chave            | QWE with Shanks transformation        | adaptive | 13                   |
    +----------------------+---------------------------------------+----------+----------------------+
    | FFTLog               | FFTLog (Hamilton 2000)                | batch    | 16                   |
    +----------------------+---------------------------------------+----------+----------------------+

.. note::

   The digital filters are named after the published filter they implement.
   The three authors are **independent designs**, so accuracy is *not* monotonic in the
   filter size: ``DHT_Key_101`` is not a refinement of ``DHT_Key_51``, and can be less
   accurate than the cheaper ``DHT_Guptasarma_Fast`` on the same integrand.
   Always compare candidates on your own form factor rather than assuming a larger
   filter is better.



.. _fixed-vs-adaptive:

Fixed, adaptive, and batch strategies
--------------------------------------

The distinction that matters most in practice is **not** the mathematical method but
whether the quadrature nodes depend on the integrand.

**Fixed-node** strategies (all the ``DHT_*`` filters and ``Fixed_DE_Ogata``) sample the
integrand at points that are known in advance:

.. math::

    G(x) = \sum_i \mathrm{base}_i \, w_i / x^2 \; f(\mathrm{base}_i / x)

They are not iterative, rely on pre-computed constants, and do not
allow the user to control the accuracy of the solution in any way - 
accuracy is set by the choice of filter.

**Adaptive** strategies (``Adaptive_DE_Ooura``, ``QWE_Key``, ``QWE_Chave``) choose their
nodes from the integrand as they go, so there is no tap matrix and the integrand must be
supplied as a callback that is invoked one point at a time. They require a tolerance
parameter (``eps_rel``): iterations stop once the difference between two consecutive
solutions falls below the tolerance. They are the right choice when the form factor is
complex or oscillatory — cases where the smooth-function assumption of FFTLog breaks down.

Note that the two double-exponential strategies fall on **opposite sides** of this
divide — ``Fixed_DE_Ogata`` places its nodes at the zeros of the Bessel function, fixed
once ``nu``, ``n_eval`` and ``f_max`` are known, whereas ``Adaptive_DE_Ooura`` is
iterative. This is why they carry an explicit ``Fixed_`` / ``Adaptive_`` qualifier while
the ``DHT_*`` (always fixed) and ``QWE_*`` (always adaptive) families do not.

**Batch (FFT-based)** strategies (``FFTLog``) are neither fixed-node nor adaptive in the
sense above. ``FFTLog`` processes the entire output array at once using a fast Fourier
transform on a log-spaced grid (Hamilton 2000, MNRAS 312, 257). Its cost is O(N log N)
across all N output points, rather than a per-point cost multiplied by N, making it
significantly faster than adaptive strategies when evaluating many points. It also
requires no ``strategy_params`` at all — no tolerance, no number of evaluations to tune.
For smooth, slowly-varying form factors on a wide dynamic range it is typically the most
accurate option as well, since adaptive quadrature accumulates numerical error point by
point whereas FFTLog processes the whole array in a single coherent pass. It is less
suited to strongly oscillatory or discontinuous form factors, and accuracy degrades at
the edges of the x array where the log-spaced grid is sparse.
