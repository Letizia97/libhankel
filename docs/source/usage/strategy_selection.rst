.. _strategy-selection:

Strategy Selection
================================

LibHankel includes various implementations of the Hankel transform, corresponding to 
different ways of approximating the Hankel transform integral. 

Independently of the interface being used, the ``hankel_transform`` function is designed 
to receive a parameter named ``strategy_name``, which allows the user the select one of 
the strategies available. 

The following table lists the strategies that are available in LibHankel.
The column named *Strategy Name* contains the possible strings the user could pass in as
``strategy_name`` when calling the ``hankel_transform`` function.
*Nodes* records whether the quadrature points are fixed or chosen adaptively (see
:ref:`fixed-vs-adaptive` below), and *n_taps* how many times the integrand is
evaluated per output point.
The column named *Index* refers to the numbering used in SASfit, for cross-reference.


.. table:: Hankel transform strategies available in LibHankel.

    +----------------------+---------------------------------------+----------+----------------+----------------------+
    | Strategy Name        | Type / Method                         | Nodes    | n_taps         | Index (as in SASfit) |
    +======================+=======================================+==========+================+======================+
    | DHT_Guptasarma       | Digital filter (Guptasarma & Singh)   | fixed    | 120 (J0)       | 6                    |
    |                      |                                       |          | 140 (J1)       |                      |
    +----------------------+---------------------------------------+----------+----------------+----------------------+
    | DHT_Guptasarma_Fast  | Digital filter (Guptasarma & Singh)   | fixed    | 61 (J0)        | 7                    |
    |                      |                                       |          | 47 (J1)        |                      |
    +----------------------+---------------------------------------+----------+----------------+----------------------+
    | DHT_Key_51           | Digital filter (Key)                  | fixed    | 51             | 8                    |
    +----------------------+---------------------------------------+----------+----------------+----------------------+
    | DHT_Key_101          | Digital filter (Key)                  | fixed    | 101            | 9                    |
    +----------------------+---------------------------------------+----------+----------------+----------------------+
    | DHT_Key_201          | Digital filter (Key)                  | fixed    | 201            | 10                   |
    +----------------------+---------------------------------------+----------+----------------+----------------------+
    | DHT_Anderson_801     | Digital filter (Anderson)             | fixed    | 801            | 11                   |
    +----------------------+---------------------------------------+----------+----------------+----------------------+
    | Fixed_DE_Ogata       | Double-exponential quadrature         | fixed    | ``n_eval``     | 1                    |
    |                      | at the zeros of the Bessel function   |          |                |                      |
    +----------------------+---------------------------------------+----------+----------------+----------------------+
    | Adaptive_DE_Ooura    | Double-exponential quadrature         | adaptive | up to          | 0                    |
    |                      |                                       |          | ``n_eval``     |                      |
    +----------------------+---------------------------------------+----------+----------------+----------------------+
    | QWE_Key              | QWE with continued fraction expansion | adaptive | up to          | 12                   |
    |                      |                                       |          | ``n_eval``     |                      |
    +----------------------+---------------------------------------+----------+----------------+----------------------+
    | QWE_Chave            | QWE with Shanks transformation        | adaptive | up to          | 13                   |
    |                      |                                       |          | ``n_eval``     |                      |
    +----------------------+---------------------------------------+----------+----------------+----------------------+

.. note::

   The digital filters are named after the published filter they implement .
   The three authors are **independent designs**, so accuracy is *not* monotonic in the
   tap count: ``DHT_Key_101`` is not a refinement of ``DHT_Key_51``, and can be less
   accurate than the cheaper ``DHT_Guptasarma_Fast`` on the same integrand.
   Always compare candidates on your own form factor rather than assuming a bigger
   filter is better.

   ``DHT_Guptasarma`` carries no tap count because theirs depends on the order of the
   Bessel function - see the *n_taps* column.

   The SASfit-indexed spellings (``DHT_6`` ... ``DHT_11``, ``DE_Ogata``, ``DE_Ooura``)
   are **no longer accepted**; use the names in the first column. The *Index* column is
   there only to let you match a strategy against its SASfit counterpart.


.. _fixed-vs-adaptive:

Fixed and adaptive strategies
-----------------------------

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
supplied as a callback that is invoked one point at a time. They require the user to
supply a tolerance parameter (``eps_rel``), used to decide when to stop iterating:
iterations stop once the difference between two consecutive solutions falls below the
tolerance. They generally perform better than the digital filters on complex or
oscillatory form factors, though they are slower.

Note that the two double-exponential strategies fall on **opposite sides** of this
divide - ``Fixed_DE_Ogata`` places its nodes at the zeros of the Bessel function, fixed
once ``nu``, ``n_eval`` and ``f_max`` are known, whereas ``Adaptive_DE_Ooura`` is
iterative. This is why they carry an explicit ``Fixed_`` / ``Adaptive_`` qualifier while
the ``DHT_*`` (always fixed) and ``QWE_*`` (always adaptive) families do not.

