.. _c-examples:

C examples
===========

Below are several example scripts on using the hankel_transform function through the :ref:`c-api`:
with a built-in form factor, with a custom input function, with a form factor supplied
as tabulated data, with the FFTLog strategy, and computing G(0) via :c:func:`hankel_G0`.

If you are unsure which strategy to start with, **FFTLog** is a good default: it requires
no strategy parameters (pass a zero-initialised ``strategy_params`` struct), runs in
O(N log N) time across the whole output array, and is typically the most accurate option
for smooth, slowly-varying form factors on a wide dynamic range. Please check the 
:ref:`c-examples-fftlog`.


To run the following examples, please install LibHankel following the instructions 
in :ref:`quickstart-c` . After that, to run the examples, you might need to refresh the 
system library cache with:

.. code-block:: bash

    sudo ldconfig 


If you also run the :ref:`python-examples`, please note that there may be slight differences
in precision with respect to the examples below, but that is expected.


.. _c-examples-builtin-f:

Example with built-in form factor
-----------------------------------

The following is an example with a built-in form factor.
Assuming you have LibHankel installed, you should be able to compile and run this example with: 

.. code-block:: bash
    
    cd libhankel
    sudo ldconfig
    gcc  examples/c/example_usage_g_dab.c -llibhankel -o example_usage_g_dab
    ./example_usage_g_dab

Here is the example itself:

.. literalinclude:: ../../../examples/c/example_usage_g_dab.c
   :language: c




.. _c-examples-custom-f:

Example with custom form factor
-----------------------------------

The following is an example with a custom form factor.
Assuming you have LibHankel installed, you should be able to compile and run this example with: 

.. code-block:: bash

    cd libhankel
    sudo ldconfig
    gcc  examples/c/example_usage_custom_form_factor.c -llibhankel -o example_usage_custom_form_factor -lm
    ./example_usage_custom_form_factor


Here is the example itself:

.. literalinclude:: ../../../examples/c/example_usage_custom_form_factor.c
   :language: c


.. _c-examples-tabulated-f:

Example with a tabulated form factor
-----------------------------------

The following example transforms a form factor that is known only at a set of
points, by handing the table to :c:func:`tabulated_ff_create`.

The interpolation itself is only half the job. Because the strategies evaluate
the form factor far outside the tabulated range, the high-:math:`q` tail has to
be defined too -- and getting that wrong fails silently rather than loudly.
This example uses a fitted power law; :ref:`tabulated-form-factors` explains
the choice.

The program is self-checking: it builds its table by sampling the built-in
``g_dab`` form factor, so it can print the transform of the tabulated data
alongside the transform of the formula it came from. The two agree to about
one part in a million.

.. code-block:: bash

    cd libhankel
    sudo ldconfig
    gcc examples/c/example_usage_tabulated_form_factor.c -llibhankel -o example_usage_tabulated_form_factor -lm
    ./example_usage_tabulated_form_factor


Here is the example itself:

.. literalinclude:: ../../../examples/c/example_usage_tabulated_form_factor.c
   :language: c


.. _c-examples-fftlog:

How to use the strategy FFTLog
-----------------------------------------

The following example runs the same ``g_dab`` form factor through the FFTLog strategy.
Unlike the adaptive strategies, FFTLog requires no convergence parameters — pass a
zero-initialised ``strategy_params`` struct and the library handles everything. It is also
faster (O(N log N) across all output points) and typically more accurate for smooth form
factors, making it the recommended starting point for most use cases.

.. code-block:: bash

    cd libhankel
    sudo ldconfig
    gcc examples/c/example_usage_fftlog.c -llibhankel -o example_usage_fftlog -lm
    ./example_usage_fftlog

Here is the example itself:

.. literalinclude:: ../../../examples/c/example_usage_fftlog.c
   :language: c


.. _c-examples-G0:

How to use hankel_G0
---------------------------------

The following example computes G(0), the Hankel transform evaluated at r = 0,
using the :c:func:`hankel_G0` function with the built-in ``g_dab`` form factor.

.. code-block:: bash

    cd libhankel
    sudo ldconfig
    gcc examples/c/example_usage_G0.c -llibhankel -o example_usage_G0 -lm
    ./example_usage_G0

Here is the example itself:

.. literalinclude:: ../../../examples/c/example_usage_G0.c
   :language: c
