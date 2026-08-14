SESANS workflow
===============

Spin-echo SANS (SESANS), as measured for example on Larmor at the ISIS Neutron
and Muon Source, does not give I(q) directly: it gives the beam polarisation as
a function of the spin-echo length :math:`\xi`. Getting from a scattering model
to something comparable with the measurement means taking a zeroth-order Hankel
transform of I(q):

.. math::

   G(\xi) = \frac{1}{2\pi} \int_0^\infty q\, J_0(q\xi)\, I(q)\, \mathrm{d}q,
   \qquad P(\xi) = G(\xi) - G(0)

This folder records the workflow that does this in SasView, so that it is clear
where LibHankel fits in, and checks that LibHankel reproduces it.

The sasmodels side
------------------

In SasView the transform is
`sasmodels.sesans.SesansTransform <https://github.com/SasView/sasmodels/blob/master/sasmodels/sesans.py>`_,
which discretises the integral above: it builds a dense matrix ``H`` of
:math:`J_0(q\xi)\, q\, \mathrm{d}q / 2\pi` on a log-spaced q grid (spacing
1.0003 by default) and applies it with a dot product.

Nothing calls ``SesansTransform`` directly. It is reached from any fit or
calculation on SESANS data through this chain:

.. code-block:: text

   DirectModel(data, model)                        direct_model.py
     DataMixin._interpret_data                     data.isSesans is True
       _make_sesans_transform(data)                direct_model.py:157
         SesansTransform(z, SElength, lam, ...)    sesans.py
           SesansTransform._set_hankel             builds H and H0

   calculator(**pars)
     DataMixin._calc_theory
       call_kernel(...)                            I(q) on transform.q_calc
       SesansTransform.apply(Iq)                   G = H.T @ Iq, G0 = H0 @ Iq
                                                   returns P = G - G0

The q range is chosen from the spin-echo lengths themselves
(:math:`q_{\min} \sim 2\pi / (N \xi_{\max})`, :math:`q_{\max} = 2\pi / \Delta\xi`),
and points where the scattering angle exceeds the acceptance of the analyser
are zeroed. There is also a one-line entry point,
``sasmodels.direct_model.Gxi(model, xi, **pars)``, which wraps the whole chain.

What the script does
--------------------

``sesans_workflow.py`` runs that workflow for a dilute dispersion of spheres of
radius 1000 A, then computes the same P(:math:`\xi`) with LibHankel doing the
transform instead of the matrix, and prints both.

The two I(q) differ by a constant factor (sasmodels normalises by particle
volume and converts to cm\ :sup:`-1`), so the comparison is made on
P(:math:`\xi`)/G(0), which does not depend on that factor.

Running it
----------

The script needs both LibHankel (see the installation instructions in the top
level ``README.rst``) and sasmodels:

.. code-block:: bash

   python -m pip install sasmodels
   python examples/sesans/sesans_workflow.py
   python examples/sesans/benchmark.py

Output of ``sesans_workflow.py``, trimmed:

.. code-block:: text

   sphere of radius 1000 A, 50 spin-echo lengths
   sasmodels : SesansTransform on 43797 q points, q = 2.51e-06 to 1.27 1/A
   libhankel : QWE_Chave over the full q range, no grid

       xi (A)      sasmodels      libhankel   difference
         50.0      -0.005476      -0.005477    -6.06e-07
        524.1      -0.273732      -0.273732    -6.27e-08
       1953.5      -0.999756      -0.999755     1.40e-06
       5000.0      -1.000001      -1.000000     1.40e-06

   largest absolute difference in P(xi)/G(0): 1.40e-06

The two agree to about 1e-6 over the whole range, and both give
P(:math:`\xi`)/G(0) = -1 beyond :math:`\xi = 2R`, where a sphere has no
remaining correlation.

The difference between the two routes is not accuracy at these settings but
what they need: sasmodels has to pick a q grid, and its result depends on that
choice, whereas LibHankel integrates over the whole half-line and takes only a
tolerance. Tested with sasmodels 1.1.0.

Speed
-----

``benchmark.py`` times the two, for 50 spin-echo lengths on a sphere of radius
1000 A. The costs have different shapes, so they are reported separately:
sasmodels pays a one-off cost to build the matrix (its size is set by the
spin-echo lengths, not by the model) and then a cost per model evaluation,
while LibHankel has no setup and does everything in the call. In a fit the
build happens once and the evaluation once per iteration, so the number to
compare against LibHankel is the per-evaluation cost.

Accuracy is in the table too, as the largest deviation of P(:math:`\xi`)/G(0)
from the analytic sphere result, because time on its own says nothing.

.. code-block:: text

   log-spaced spin-echo lengths, n = 50
     sasmodels  build SesansTransform    45.93 ms  (43797 q points)
     sasmodels  per model evaluation      7.36 ms  (best 2.87)   max error 1.4e-06
     libhankel strategy        median      best  vs sasmodels   max error
     DHT_Guptasarma_Fast      0.08 ms   0.07 ms         0.01x     5.8e-03
     DHT_Key_51               0.06 ms   0.06 ms         0.01x     1.1e-02
     DHT_Key_201              0.22 ms   0.21 ms         0.03x     5.8e-04
     DHT_Anderson_801         1.45 ms   1.44 ms         0.20x     2.0e-04
     Fixed_DE_Ogata           0.52 ms   0.51 ms         0.07x     5.3e-04
     Adaptive_DE_Ooura        5.50 ms   3.64 ms         0.75x     8.2e-06
     QWE_Key                  5.65 ms   5.58 ms         0.77x     1.4e-07
     QWE_Chave                6.01 ms   5.87 ms         0.82x     1.4e-07

   linear spin-echo lengths, n = 50
     sasmodels  build SesansTransform    21.83 ms  (33727 q points)
     sasmodels  per model evaluation      3.98 ms  (best 2.09)   max error 3.1e-04
     libhankel strategy        median      best  vs sasmodels   max error
     DHT_Guptasarma_Fast      0.05 ms   0.05 ms         0.01x     3.2e-03
     DHT_Key_51               0.06 ms   0.05 ms         0.01x     1.1e-02
     DHT_Key_201              0.20 ms   0.20 ms         0.05x     2.9e-04
     DHT_Anderson_801         1.41 ms   1.40 ms         0.36x     2.2e-04
     Fixed_DE_Ogata           0.50 ms   0.49 ms         0.13x     5.3e-04
     Adaptive_DE_Ooura        4.83 ms   2.65 ms         1.21x     1.1e-05
     QWE_Key                  9.13 ms   9.00 ms         2.30x     8.9e-08
     QWE_Chave               19.43 ms  18.24 ms         4.88x     8.9e-08

Reading that:

* The digital filters are 30 to 100 times faster than a sasmodels evaluation,
  at 1e-2 to 1e-3 accuracy. ``DHT_Key_201`` is the useful point on that curve:
  0.2 ms at 3e-4, which matches what sasmodels itself achieves on the linear
  grid, about 20 times faster.
* The adaptive strategies land in the same range as sasmodels per evaluation
  (0.6x to 5x, depending on grid and tolerance) and buy three to four more
  digits. ``QWE_Chave`` is the slowest and the most accurate.
* The sasmodels evaluation is memory-bound rather than compute-bound: the H
  matrix here is 43797 x 50 doubles, about 17 MB read per evaluation. That is
  why its timings scatter (median 7.4 ms, best 2.9 ms) while the LibHankel
  timings barely move, and why its cost grows with the density of the
  spin-echo grid rather than with the complexity of the model.
* Nothing here is amortised on the LibHankel side. Its cost is per call and
  scales with the number of spin-echo lengths, whereas sasmodels' 22 to 46 ms
  build is paid once per dataset.

Two things to be aware of when reproducing this:

* At the default ``n_eval`` of 250, ``QWE_Key`` and ``QWE_Chave`` fail to
  converge on the linear grid for the two :math:`\xi` values either side of
  2R, where G(:math:`\xi`) passes through zero and a relative tolerance cannot
  be met. ``benchmark.py`` uses ``n_eval = 500``, ``eps_rel = 1e-8``.
* ``Fixed_DE_Ogata`` needs ``f_max`` around 1e-4 for this integrand; at
  ``f_max = 1`` it returns an answer that is wrong by 100%, without failing.

Timings were taken on an Intel Core Ultra 7 165H under WSL2, best and median of
15 calls. Treat the ratios as indicative.

Speed on real data
------------------

``benchmark.py`` invents its spin-echo lengths, and the shape of that grid is
exactly what sets the cost of the sasmodels route.
``benchmark_real_data.py`` repeats the comparison on measured SESANS files
instead. It uses the three that ship with ``sasdata``, so there is nothing to
download:

.. code-block:: bash

   python -m pip install sasmodels sasdata
   python examples/sesans/benchmark_real_data.py

What the data supplies is the spin-echo lengths, the wavelengths and the
analyser acceptance, which between them fix the q grid and the mask. It does
not supply I(q): the transform is always applied to model I(q) on the grid the
transform chose, and the measured polarisation is the fit target, not an input.
So these runs read real files, build the transform they imply, and put a sphere
model through it.

.. code-block:: text

   sphere_isis.ses   PMMA in deuterated decalin, ISIS, time-of-flight
     57 spin-echo lengths, 260 to 19303 A;  lambda 1.61 to 13.89 A;  theta_max 0.09 rad
     sasmodels  build SesansTransform    52.38 ms  (43936 q points, 20.0 MB, 13.1% of H zeroed)
     sasmodels  per model evaluation      8.33 ms  (best 4.37)   max error 1.0e-05
     libhankel strategy        median      best  vs sasmodels   max error
     DHT_Key_51               0.07 ms   0.07 ms         0.01x     1.5e-03
     DHT_Key_201              0.28 ms   0.25 ms         0.03x     2.8e-04
     DHT_Anderson_801         1.71 ms   1.59 ms         0.20x     2.4e-04
     Fixed_DE_Ogata           0.60 ms   0.53 ms         0.07x     3.5e-03
     Adaptive_DE_Ooura        3.81 ms   2.41 ms         0.46x     8.1e-06
     QWE_Key                  7.63 ms   7.18 ms         0.92x     5.7e-07
     QWE_Chave               14.09 ms  13.00 ms         1.69x     5.7e-07

   sphere2micron.ses   2 um polystyrene, Delft, monochromatic
     40 spin-echo lengths, 392 to 46099 A;  lambda 2.11 A;  theta_max 0.0168 rad
     sasmodels  build SesansTransform    17.96 ms  (32216 q points, 10.3 MB, 0.0% of H zeroed)
     sasmodels  per model evaluation      5.07 ms  (best 3.95)   max error 3.8e-04
     libhankel strategy        median      best  vs sasmodels   max error
     DHT_Key_51               0.05 ms   0.05 ms         0.01x     1.3e-02
     DHT_Key_101              0.09 ms   0.09 ms         0.02x     2.7e-02
     DHT_Key_201              0.17 ms   0.16 ms         0.03x     3.5e-04
     DHT_Anderson_801         1.28 ms   1.19 ms         0.25x     3.3e-04
     Fixed_DE_Ogata           0.48 ms   0.47 ms         0.10x     3.1e-04
     Adaptive_DE_Ooura        4.17 ms   2.57 ms         0.82x     3.3e-06
     QWE_Key                        -         -             -   no result: 1/40 points fail
     QWE_Chave                      -         -             -   no result: 1/40 points fail

   sphere2micron_long.ses   the same sample, out to 20 um
     84 spin-echo lengths, 392 to 200803 A;  lambda 2.11 A;  theta_max 0.0168 rad
     sasmodels  build SesansTransform    54.59 ms  (39595 q points, 26.6 MB, 0.0% of H zeroed)
     sasmodels  per model evaluation      6.73 ms  (best 3.79)   max error 3.8e-04
     libhankel strategy        median      best  vs sasmodels   max error
     DHT_Key_201              0.34 ms   0.33 ms         0.05x     3.5e-04
     DHT_Anderson_801         2.48 ms   2.42 ms         0.37x     3.3e-04
     Fixed_DE_Ogata           0.68 ms   0.65 ms         0.10x     3.5e-03
     Adaptive_DE_Ooura        3.02 ms   2.74 ms         0.45x     3.3e-06
     QWE_Key                        -         -             -   no result: 1/84 points fail

The picture is the same as on the synthetic grids, with three things that only
show up on real data:

* **The QWE strategies return nothing at all on two of the three datasets.**
  ``sphere2micron.ses`` has a spin-echo length at 20316 A, just past
  2R = 20000 A, where G(:math:`\xi`) crosses zero and a relative tolerance
  cannot be met. That single point fails, and because ``hankel_transform``
  raises for the whole array, the other 39 points are lost with it. This is the
  same failure the synthetic benchmark hits, but there it is a curiosity and
  here it takes out the most accurate strategy on a real measurement. A
  per-point fallback, or an absolute-tolerance floor, would be needed to use
  QWE on data like this.
* **The digital filters are not monotonic in filter length.** On the 2 um data
  ``DHT_Key_101`` (2.7e-02) is worse than ``DHT_Key_51`` (1.3e-02), while on the
  ISIS data the ordering is the expected one. Filter length is not a dial for
  accuracy on a given integrand.
* With QWE out, ``Adaptive_DE_Ooura`` is the most accurate strategy that
  survives every real grid, and it is still 1.2 to 2.2 times faster than a
  sasmodels evaluation. Its error is at or below what sasmodels itself achieves
  (8.1e-06 against 1.0e-05 on the ISIS grid, 3.3e-06 against 3.8e-04 on the
  2 um grids). ``DHT_Key_201`` matches the sasmodels error exactly on those
  grids, at 0.17 to 0.34 ms, some 20 to 30 times faster than an evaluation.

The ISIS file is the only one where the analyser acceptance does anything:
13.1% of the H matrix is zeroed, because it is a time-of-flight measurement and
the long-wavelength frames scatter past the analyser. Even so, sasmodels lands
1.0e-05 from the ideal-transform answer and LibHankel 5.7e-07, both far inside
the measurement noise. Worth noting that the mask is a physical correction
LibHankel has no equivalent of: on this dataset it is negligible, but that is a
property of the dataset, not a general result.

Parameters are nominal (the quoted particle size, unit contrast, no
background), since the comparison is on P(:math:`\xi`)/G(0) and neither
contrast nor scale survives that normalisation. Only the radius matters, and
these are calibration samples chosen precisely because it is known.

Speed on tabulated I(q)
-----------------------

Everything above transforms an analytic form factor, which either route can
evaluate wherever it likes. ``benchmark_tabulated.py`` transforms a *table* of
(q, I) pairs instead, which is the shape scattering data actually arrives in,
and it reverses most of the conclusions.

The two routes need a table in different ways:

* ``SesansTransform`` never picks a q of its own. ``apply(Iq)`` is
  ``H.T @ Iq``, so it needs I(q) once, on the fixed grid ``q_calc``. Getting a
  table onto that grid is one vectorised interpolation and the transform is
  otherwise untouched.
* LibHankel picks its own abscissae, so the table has to be wrapped in a
  callable. Through the Python interface that is one ``PyObject_CallObject``
  per evaluation (``src/py_interface.c:20-36``, which also rebuilds the
  parameter list on every call), and the strategy decides how many.

.. code-block:: text

   sphere radius 1000 A, 57 spin-echo lengths from sphere_isis.ses
   tabulated I(q): 1000 points, 1e-06 to 10 1/A, log-log interpolated
     sasmodels  interpolate + apply      10.30 ms  (best 4.76)   max error 4.1e-04
     libhankel strategy        median      best  vs sasmodels   max error   I(q) calls   outside
     DHT_Guptasarma_Fast      3.73 ms   3.08 ms         0.36x     3.0e-03         3477       35%
     DHT_Key_51               5.39 ms   4.73 ms         0.52x     1.6e-03         2907        2%
     DHT_Key_201              4.77 ms   4.69 ms         0.46x     3.4e-04        11457       35%
     DHT_Anderson_801        11.45 ms  11.10 ms         1.11x     2.5e-04        45657       80%
     Fixed_DE_Ogata          10.50 ms   8.91 ms         1.02x     3.5e-03        14250        1%
     Adaptive_DE_Ooura       67.85 ms  62.27 ms         6.58x     1.6e-04       108289        2%
     QWE_Key                        -         -             -   no result       508050        0%
     QWE_Chave                      -         -             -   no result       539501        0%

Against the same strategies on the same :math:`\xi` grid with the built-in C
sphere, this is a 7 to 18 times slowdown: ``DHT_Key_201`` goes from 0.28 ms to
4.77 ms, ``Adaptive_DE_Ooura`` from 3.81 ms to 67.85 ms. The 30-to-100-times
advantage over sasmodels becomes at best about two times, and the sasmodels
route is now the faster one for half the strategies.

* **The QWE strategies fail outright**, after half a million interpolator calls
  and about a second of work. Log-log interpolation makes the integrand
  piecewise smooth with a kink at every tabulated point, and the adaptive
  bisection keeps subdividing against those kinks. This is a different failure
  from the zero-crossing one above, and it rules the QWE strategies out for
  tabulated input entirely rather than for particular :math:`\xi`.
* **The digital filters have an accuracy floor that more data does not lift.**
  Refining the table stops helping them at about 1000 points, while sasmodels
  keeps improving:

  .. code-block:: text

       points    sasmodels   DHT_Key_201   DHT_Anderson_801
          100      1.7e-02       1.1e-02            1.2e-02
          300      3.0e-03       1.7e-03            1.8e-03
         1000      4.1e-04       3.4e-04            2.5e-04
         3000      7.0e-05       2.9e-04            2.4e-04
        10000      1.6e-05       3.0e-04            2.4e-04

  Below about 1000 points the filters are slightly ahead; above it they
  saturate at their own quadrature error while the dense-matrix route converges
  on the table. For a well-sampled curve sasmodels ends up an order of
  magnitude more accurate.
* **The wide filters waste most of their samples.** ``DHT_Anderson_801`` spends
  80% of its 45657 evaluations outside the table, where the answer is zero by
  definition. Its nodes scale as 1/:math:`\xi` over ten decades, far wider than
  any table covers, so it pays full price for samples that cannot contribute.

None of that is the main limitation, though. The q window is:

.. code-block:: text

   window                          q range (1/A)    sasmodels   DHT_Key_201
   D22 as shipped           4.0e-02 to  6.1e-01      9.2e-01       1.0e+00
   typical SANS             3.0e-03 to  5.0e-01      6.4e-01       6.7e-01
   SANS + USANS             1.0e-04 to  5.0e-01      2.2e-03       2.6e-03
   USANS to low q           1.0e-05 to  5.0e-01      2.1e-04       2.6e-04
   generous synthetic       1.0e-06 to  1.0e+01      4.1e-04       3.4e-04

A 1000 A sphere seen through a real SANS window gives an answer that is wrong
by 100%, on both routes equally, because the q that carries the SESANS signal
is simply not measured: ``q_calc`` for this dataset runs from 5.7e-07 to
3.0e-01 1/A, and D22 starts at 4.0e-02. Useful accuracy needs the table to
reach 1e-4, and matching the model-based results needs 1e-5 — USANS territory.
This is a property of the measurement, not of either transform, and it is why
SESANS analysis fits a model rather than transforming a measured curve.

So the honest summary is that the choice of transform barely matters on
tabulated data. If it has to be made, ``DHT_Key_201`` is again the sensible
one: about twice as fast as the sasmodels route, comparable error at realistic
table densities, and no failure modes. The larger point is that the C-level
form factor is where LibHankel's advantage lives, and a Python callback gives
most of it back — a batched interface, handing the strategy's abscissae out and
taking a vector of I(q) values in, would recover it, and the fixed-node filters
are the strategies that could offer one.

Replacing the transform
=======================

That last observation is worth following, because sasmodels already works that
way. Its transform never calls the model: it publishes ``q_calc``, someone else
evaluates, and ``apply`` consumes the vector. So a strategy that can say in
advance where it needs I(q) does not need a callback at all — it can *be* a
sasmodels resolution object, and the border between C and Python is crossed
once at setup instead of once per point.

Only the fixed-abscissa filters can do this. Their abscissae are a table
divided by :math:`\xi` (``src/hankel_DHT.c:110-113``), known before ``f`` is
evaluated; the adaptive strategies choose theirs from values they have already
seen, so they cannot answer the question ``q_calc`` asks. This is the same
property the batched interface would have relied on, used from the other side.

``dht_nodes`` publishes the tables:

.. code-block:: python

   abscissae, weights = libhankel.dht_nodes("DHT_Key_201", 0)

``dht_sesans.py`` builds a transform on them, and
``benchmark_dht_transform.py`` puts it through ``DirectModel`` on the same
three files. The swap is one call, and nothing downstream changes:

.. code-block:: python

   use_dht_transform("DHT_Key_201")
   calculator = DirectModel(data, model)   # now uses the filter
   calculator(radius=1000.0, **pars)       # returns P(xi) as before

.. code-block:: text

   sphere_isis.ses   57 spin-echo lengths, 260 to 19303 A;  radius 1000 A
     transform                build  per eval      best  q points  speedup   max error
     sasmodels dense       54.91 ms   6.39 ms   3.00 ms     43936     1.0x     1.0e-05
     DHT_Key_51             0.23 ms   0.12 ms   0.10 ms      2907    53.8x     1.5e-03
     DHT_Key_101            0.27 ms   0.15 ms   0.14 ms      5757    41.6x     7.4e-04
     DHT_Key_201            0.46 ms   0.25 ms   0.23 ms     11457    25.2x     2.8e-04
     DHT_Anderson_801       2.03 ms   1.42 ms   1.37 ms     45657     4.5x     2.4e-04

The other two datasets behave the same way: ``DHT_Key_201`` runs 30 to 59 times
faster per evaluation at 3.5e-04 against the dense route's 3.8e-04, and builds
in about 0.55 ms against 22 to 57 ms. The speedup figures move by a third
between runs, because the baseline is memory-bound and scatters — its median is
roughly twice its best — so read them as "tens of times", not to two figures.
The error figures are stable, and they agree with the ones
``benchmark_real_data.py`` gets by calling ``hankel_transform`` directly, which
is a useful check that the reimplementation in numpy is faithful.

Two things had to be settled to make this a real replacement rather than a
demonstration.

**G(0).** ``apply`` returns :math:`G(\xi) - G(0)`, and G(0) is the transform at
:math:`\xi = 0`. A filter evaluates at ``node/xi``, so it cannot go there;
``SesansTransform`` gets G(0) free from its dense grid, and
``benchmark_real_data.py`` sidesteps the problem by using the closed form for a
sphere, which does not exist for a general model. But G(0) is a plain integral,
not a Hankel transform, and the union of the abscissae over all the spin-echo
lengths already covers q densely — 11457 points across 13 decades — because
each :math:`\xi` shifts the same table by a different :math:`1/\xi`. A
trapezoid rule in log q over that union costs no extra model evaluations, and
it is more accurate than the rectangle rule it replaces:

.. code-block:: text

   route                 q points  relative error
   reference               400000   400k point log grid
   sasmodels dense          43936        -1.6e-04
   DHT_Key_51                2907        -1.1e-05
   DHT_Key_101               5757        -4.4e-07
   DHT_Key_201              11457         1.1e-07
   DHT_Anderson_801         45657        -4.7e-07

The one case it does not cover is a single spin-echo length, where the union is
just the table and the sampling drops to about 19 points per decade.

**The acceptance mask.** ``benchmark_real_data.py`` never applies it, so on the
ISIS file it and sasmodels are integrating slightly different things.
``DHTSesansTransform`` applies the same masking as
``SesansTransform._set_hankel``, elementwise on its own abscissae, which is
cheaper than masking a matrix because the abscissae vary with :math:`\xi`
anyway.

One incidental finding, worth knowing before copying ``apply`` from sasmodels:
writing the G(0) reduction as ``np.dot`` costs 2.9 ms on 11457 points against
0.015 ms for an explicit multiply and sum. ``np.dot`` goes to threaded BLAS,
and at this size the thread setup dwarfs the arithmetic — ``np.dot`` of two
vectors of ones is just as slow, so it is overhead, not the values. It is right
for the ``(n_q, n_xi)`` matrix product sasmodels uses it for and wrong for a
vector reduction, and the first version of this transform was 20 times slower
than it needed to be for exactly that reason.

What this does not buy is accuracy. On the ISIS grid the dense route is still
the better one, 1.0e-05 against 2.8e-04, and the shorter filters are erratic:
``DHT_Key_101`` is worse than ``DHT_Key_51`` on both 2 :math:`\mu`\ m datasets
(2.7e-02 against 1.3e-02), so filter length is not a dial that can be turned
for accuracy. ``DHT_Key_201`` is the one to use — it is the shortest filter
that is accurate on all three datasets, and it is the point where the remaining
cost has moved from the transform into the model evaluation itself.
