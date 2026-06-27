.. libhankel documentation master file, created by
   sphinx-quickstart on Wed Jun 17 17:22:03 2026.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

libhankel documentation
=======================

.. Add your content using ``reStructuredText`` syntax. See the
.. `reStructuredText <https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html>`_
.. documentation for details.

Welcome to the documentation for **Libhankel**.

LibHankel is a C library with a number of methods for computing Hankel transforms, 
re-implemented based on the strategies available in [SASfit](https://github.com/SASfit/SASfit). 
The methods are effectively different ways of approximating the integral that constitutes 
the Hankel transform.


Overview
---------------

.. toctree::
   :maxdepth: 2
   :caption: Overview

   overview/the_hankel_transform


Getting Started
---------------

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   getting_started/quickstart_python
   getting_started/quickstart_c





Examples
--------

Browse practical examples demonstrating common workflows.

.. toctree::
   :maxdepth: 1
   :caption: Examples

   examples/python_examples
   examples/c_examples



Usage
-----------

Usage info.

.. toctree::
   :maxdepth: 2
   :caption: Usage

   guides/usage



API Reference
-------------

Detailed reference documentation for both interfaces.

.. toctree::
   :maxdepth: 2
   :caption: API References

   api/python_api
   api/c_api

