.. _c-examples:

C examples
===========

Below are two example scripts on using the hankel_transform function through the :ref:`c-api`,
both with a built-in and a custom input function.  Please note that both examples give 
the same (or very similar) result. The main difference between the two is the fact 
that one uses a built-in function called "gdab", while the other defines a custom C 
function essentially containing the same code as the built-in "gdab" function.

To run the following examples, please install LibHankel following the instructions 
in :ref:`quickstart-c` . After that, to run the examples, you might need to refresh the 
system library cache with:

.. code-block:: bash

    sudo ldconfig 



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
   :lines: 1-1000




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
   :lines: 1-1000
