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
    