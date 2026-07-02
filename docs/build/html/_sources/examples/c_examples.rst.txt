.. _c-examples:

C examples
===========

Below are two example scripts on using the hankel_transform function using the C API,
both with a built-in and a custom input function.

To run the following examples, please install LibHankel following the instructions 
in :ref:`quickstart-c` . After that, you might need to refresh the 
system library cache with:

.. code-block:: bash

    sudo ldconfig 



.. _c-examples-builtin-f:

With built-in form factor
-----------------------------------

The following is an example with a built-in form factor.

.. literalinclude:: ../../../examples/c/example_usage_g_dab.c
   :language: c
   :lines: 1-1000

Assuming you have installed libhankel as mentioned at the start of this page,
you should be able to compile and run the above example with: 

.. code-block:: bash
    
    cd libhankel
    gcc  examples/c/example_usage_g_dab.c -llibhankel -o example_usage_g_dab
    ./example_usage_g_dab



.. _c-examples-custom-f:

With custom form factor
-----------------------------------

The following is an example with a custom form factor.

.. literalinclude:: ../../../examples/c/example_usage_custom_form_factor.c
   :language: c
   :lines: 1-1000

Assuming you have installed libhankel as mentioned at the start of this page,
you should be able to compile and run the above example with: 

.. code-block:: bash

    cd libhankel
    gcc  examples/c/example_usage_custom_form_factor.c -llibhankel -o example_usage_custom_form_factor -lm
    ./example_usage_custom_form_factor
