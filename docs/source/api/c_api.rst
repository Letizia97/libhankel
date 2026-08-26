.. _c-api:

C API
===========

This page contains explanations of the main functions / structs
a user might need, when using the library from C.

Please refer to the :ref:`c-examples` for details on how to use the 
hankel transform with either a built-in form factor or a custom input function.

.. _strategy_params_c_api:

Strategy params
------------------
.. doxygenstruct:: strategy_params
   :members:


.. _hankel_transform_c_api:

Hankel transform function
---------------------------

.. doxygenfunction:: hankel_transform


.. _form_factors_c_api:

For form factors
----------------------
.. doxygentypedef:: form_factor_f
.. doxygenstruct:: form_factor_ctx
   :members:
.. doxygenfunction:: form_factor_g_dab
.. doxygenfunction:: form_factor_sphere
.. doxygenfunction:: form_factor_broad_peak


.. _interpolation_c_api:

Interpolation
----------------------

For transforming a form factor that is known only at a set of points rather
than as a formula. These functions build a PCHIP cubic spline through the
data; the spline is monotonicity-preserving, so it will not overshoot into
negative intensities the way a classical C2 spline can.

.. warning::

    :c:func:`cubic_interp_eval` returns NaN outside the tabulated range, and
    the strategies evaluate the form factor far outside it. Passing the bare
    interpolator to :ref:`hankel_transform <hankel_transform_c_api>` fails
    silently -- with a status code of 0. Use
    :ref:`tabulated_ff <tabulated_ff_c_api>` below, which handles the range
    outside the table for you, unless you have a reason to interpolate
    directly.

.. doxygentypedef:: cubic_interp_t
.. doxygenfunction:: cubic_interp_create
.. doxygenfunction:: cubic_interp_eval
.. doxygenfunction:: cubic_interp_destroy


.. _tabulated_ff_c_api:

Tabulated form factors
----------------------------

A form factor defined by a table of points rather than by a formula, ready to
hand to :ref:`hankel_transform <hankel_transform_c_api>`. This is the
interpolation above plus the rules for what happens outside the table, which
is the part that has to be right: the strategies sample over tens of decades
in :math:`q`, so every call goes outside the tabulated range.

:c:func:`tabulated_ff_eval` has the signature of a :c:type:`form_factor_f`,
and a handle is its own context, so it is passed straight through as the
``f_ctx`` argument with no adapter to write:

.. code-block:: c

    tabulated_ff_t *ff = NULL;
    int status = tabulated_ff_create(q, f, n, TABULATED_TAIL_POWER_LAW, 0.0, &ff);
    if (status != 0) { /* see the status codes page */ }

    strategy_params params = {.n_eval = 250, .eps_rel = 1e-9};
    status = hankel_transform(0, tabulated_ff_eval, x, n_x, ff, result, "QWE_Key", params);

    tabulated_ff_destroy(ff);

Passing ``0.0`` as the exponent fits it from the table. See
:ref:`tabulated-form-factors` for how to choose the tail, and
:ref:`c-examples-tabulated-f` for a complete program.

.. doxygenenum:: tabulated_tail_t
.. doxygentypedef:: tabulated_ff_t
.. doxygenfunction:: tabulated_ff_create
.. doxygenfunction:: tabulated_ff_eval
.. doxygenfunction:: tabulated_ff_destroy
