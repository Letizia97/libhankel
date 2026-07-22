
.. _build-docs:

Building the docs locally
==========================

Official docs are available 
`through github.io <https://letizia97.github.io/libhankel/index.html>`_ . 
These are based on the "main" branch. When working on a dev branch, then 
docs must be built locally to check changes. Once the dev branch is merged into
main, the github.io docs will be re-deployed automatically to reflect the merged changes.

To build the docs locally:

.. code-block:: bash

   cd libhankel
   pip install -e .
   doxygen Doxyfile
   cd libhankel/docs
   make html

