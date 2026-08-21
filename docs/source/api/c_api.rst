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
    silently -- with a status code of 0. Read
    :ref:`tabulated-form-factors` before using these.

.. doxygentypedef:: cubic_interp_t
.. doxygenfunction:: cubic_interp_create
.. doxygenfunction:: cubic_interp_eval
.. doxygenfunction:: cubic_interp_destroy
    