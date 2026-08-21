.. _c-examples:

C examples
===========

Below are three example scripts on using the hankel_transform function through the :ref:`c-api`:
with a built-in form factor, with a custom input function, and with a form factor supplied
as tabulated data. Please note that all three examples give the same (or very similar) result.
The first uses a built-in function called "gdab"; the second defines a custom C function
essentially containing the same code as the built-in "gdab" function; the third interpolates
a table of values sampled from it.

To run the following examples, please install LibHankel following the instructions 
in :ref:`quickstart-c` . After that, to run the examples, you might need to refresh the 
system library cache with:

.. code-block:: bash

    sudo ldconfig 


If you also run the :ref:`python-examples`, please note that there may be slight differences
in precision with respect to the examples below, but that is expected.


.. _c-examples-builtin-f:

With built-in form factor
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

With custom form factor
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

With a tabulated form factor
-----------------------------------

The following example transforms a form factor that is known only at a set of
points, by building a PCHIP spline through it with the
:ref:`interpolation functions <interpolation_c_api>`.

The interpolation itself is only half the job. Because the strategies evaluate
the form factor far outside the tabulated range, the callback also has to
define the two tails -- and getting that wrong fails silently rather than
loudly. :ref:`tabulated-form-factors` explains the reasoning; this example
implements it.

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
