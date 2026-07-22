.. libhankel documentation master file, created by
   sphinx-quickstart on Wed Jun 17 17:22:03 2026.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

LibHankel documentation
=======================

.. Add your content using ``reStructuredText`` syntax. See the
.. `reStructuredText <https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html>`_
.. documentation for details.

Welcome to the documentation for **Libhankel**.

LibHankel is a C library with a number of methods for computing Hankel transforms, 
re-implemented based on the strategies available in `SASfit <https://github.com/SASfit/SASfit>`_. 
The methods are effectively different ways of approximating the integral that constitutes 
the Hankel transform.



.. toctree::
   :maxdepth: 2
   :titlesonly:
   :caption: Overview

   overview/the_hankel_transform



.. toctree::
   :maxdepth: 2
   :titlesonly:
   :caption: Getting Started

   getting_started/quickstart_python
   getting_started/quickstart_c



.. toctree::
   :maxdepth: 1
   :caption: Examples

   examples/python_examples
   examples/c_examples



.. toctree::
   :maxdepth: 2
   :caption: Usage

   usage/strategy_selection
   usage/strategy_params
   usage/status_codes



.. toctree::
   :maxdepth: 2
   :caption: API References
   
   api/python_api
   api/c_api


.. toctree::
   :maxdepth: 2
   :caption: Contributors
   
   contributors/build_docs
   contributors/installation_for_development
