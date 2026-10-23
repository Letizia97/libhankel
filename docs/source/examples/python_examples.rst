
.. _python-examples:

Python examples
================


Below are several examples on how to use the hankel_transform function
through the :ref:`python-api`, covering different ways to supply a form factor and
different strategy choices.

If you are unsure which strategy to start with, **FFTLog** is a good default: it requires
no strategy parameters (just pass an empty dict), runs in O(N log N) time across the whole
output array, and is typically the most accurate option for smooth, slowly-varying form
factors on a wide dynamic range. Please check :ref:`python-examples-fftlog` for an 
example that uses FFTLog.
The adaptive strategies (``QWE_*``, ``Adaptive_DE_Ooura``)
are better when the form factor is complex or oscillatory and you need fine control over
the convergence tolerance.

In order to run these examples, please install LibHankel through
the instructions at :ref:`quickstart-python`.

If you also run the :ref:`c-examples`, please note that there may be slight differences
in precision with respect to the examples below, but that is expected.


.. _python-examples-builtin-f:

Example with built-in form factor
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

Example with custom form factor
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


.. _python-examples-tabulated-ff:

Example with tabulated form factor
------------------------------------

The following is an example with a tabulated form factor, where the form factor
values are provided as arrays of q and f values rather than as an analytic function.

Assuming you have LibHankel installed, you should be able to run this example with:

.. code-block:: bash

   cd libhankel
   python examples/python/example_usage_tabulated_ff.py


Here is the content of the example:

.. literalinclude:: ../../../examples/python/example_usage_tabulated_ff.py
   :language: python
   :lines: 1-1000



.. _python-examples-fftlog:

Example on how to use the strategy FFTLog
-----------------------------------------

The following example runs the built-in ``gdab`` form factor through the FFTLog strategy.
Unlike the adaptive strategies, FFTLog requires no convergence parameters — pass an empty
dict ``{}`` and the library handles everything. It is also faster (O(N log N) across all
output points rather than per-point quadrature) and typically more accurate for smooth form
factors, making it the recommended starting point for most use cases.

Assuming you have LibHankel installed, you should be able to run this example with:

.. code-block:: bash

   cd libhankel
   python examples/python/example_usage_fftlog.py


Here is the content of the example:

.. literalinclude:: ../../../examples/python/example_usage_fftlog.py
   :language: python
   :lines: 1-1000


.. _python-examples-G0:

Example on how to use hankel_G0
-------------------------------

The following example computes G(0), the Hankel transform evaluated at r = 0,
using the hankel_G0 function.

Assuming you have LibHankel installed, you should be able to run this example with:

.. code-block:: bash

   cd libhankel
   python examples/python/example_usage_G0.py


Here is the content of the example:

.. literalinclude:: ../../../examples/python/example_usage_G0.py
   :language: python
   :lines: 1-1000