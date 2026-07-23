
.. _python-examples:

Python examples
================


.. _python-examples-builtin-f:


Below are a couple of examples on how to use the hankel_transform function 
through the :ref:`python-api`, both with a built-in and a custom input function. 
Please note that both examples give the same (or very similar) result. The
main difference between the two is the fact that one uses a built-in function called
"gdab", while the other defines a custom Python function essentially containing the same
code as the built-in "gdab" function.

In order to run these examples, please install LibHankel through
the instructions at :ref:`quickstart-python`. 


With built-in form factor
--------------------------

The following is an example with a built-in form factor.
Assuming you have LibHankel installed, you should be able to compile and run this example with: 

.. code-block:: bash
   
   cd libhankel
   python examples/python/example_usage_g_dab.py


Here is the content of the example:

.. literalinclude:: ../../../examples/python/example_usage_g_dab.py
   :language: python
   :lines: 1-1000




.. _python-examples-custom-f:

With custom form factor
------------------------

The following is an example with a custom form factor.
Assuming you have LibHankel installed, you should be able to compile and run this example with: 

.. code-block:: bash
   
   cd libhankel
   python examples/python/example_usage_custom_form_factor.py


Here is the content of the example:


.. literalinclude:: ../../../examples/python/example_usage_custom_form_factor.py
   :language: python
   :lines: 1-1000

