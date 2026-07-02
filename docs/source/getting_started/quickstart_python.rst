.. _quickstart-python:

Installation for Python users
=================================

The easiest way to install LibHankel is by using the Python package manager, `pip <https://pip.pypa.io/en/stable/>`__.
Alternatively, it can be installed from source through git clone. In either case, it is recommended to create a virtual
environment in which to install LibHankel, by typing for example :

.. code-block:: bash
   python3 -mvenv venv
   . venv/bin/activate


Installing via ``pip``
----------------------

LibHankel can be installed via the command line by entering:

.. code-block:: bash

   python -m pip install libhankel


Installing from source
----------------------

You may instead wish to install from source, e.g., to get the very latest version
of the code that is still in development.

1. Download this repository or clone it using:
   ``git clone https://github.com/Letizia97/libhankel.git``

2. Open up a terminal (command prompt) and go into the
   ``libhankel`` directory.

3. Once you are in the right directory, type

   .. code-block:: bash

      python -m pip install .
