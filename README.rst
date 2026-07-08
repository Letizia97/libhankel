This library is based on SASfit (https://github.com/SASfit/SASfit), whose main contributors are Joachim Kohlbrecher and Ingo Breßler.

The idea with the library is to re-implement the Hankel transform strategies available in SASfit, and provide APIs so they can be called from codebases in different languages, e.g. Python, as well as C (and therefore by SASfit itself).

Installation for C users
=================================

.. c-installation-start

LibHankel has been tested on Linux, although it should work on other platforms. It requires:

- Meson

- A C compiler


To build and install:

.. code-block:: bash

   git clone https://github.com/libhankel
   cd libhankel
   meson setup build
   meson compile -C build
   meson install -C build


.. c-installation-end



Installation for Python users
=================================

.. python-installation-start

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


.. python-installation-end

