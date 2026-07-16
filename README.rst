This library is based on SASfit (https://github.com/SASfit/SASfit), whose main contributors are Joachim Kohlbrecher and Ingo Breßler.

The idea with the library is to re-implement the Hankel transform strategies available in SASfit, and provide APIs so they can be called from codebases in different languages, e.g. Python, as well as C (and therefore by SASfit itself).

Documentation
===============

Please note, docs are available at https://letizia97.github.io/libhankel/index.html

Installation for C users
=================================

.. c-installation-start

LibHankel has been tested on Linux, although it should work on other platforms. It requires:

- Meson >= 1.4.0
- Ninja
- A C compiler (e.g. gcc)

It also requires Boost development package as a dependency. This can be installed with:

.. code-block:: bash

   sudo apt update
   sudo apt install libboost-all-dev

To build and install LibHankel, use:

.. code-block:: bash

   git clone https://github.com/Letizia97/libhankel.git
   cd libhankel
   meson setup build
   meson compile -C build
   meson install -C build


.. c-installation-end

You can check the installation by running the C examples available through the docs.



Installation for Python users
=================================

.. python-installation-start

The easiest way to install LibHankel is by using the Python package manager, `pip <https://pip.pypa.io/en/stable/>`__.
Alternatively, it can be installed from source through git clone.


Installing via ``pip``
----------------------

**INSTALLATION THROUGH PIP IS NOT YET AVAILABLE, PLEASE IGNORE THIS SECTION AND INSTALL FROM SOURCE 
(FOLLOW INSTRUCTIONS IN SECTION "Installing from source" BELOW).**

It is recommended to first create a virtual environment in which to install LibHankel, e.g., 
through the command line:

.. code-block:: bash
   
   python3 -mvenv venv
   . venv/bin/activate


LibHankel can then be installed by entering:

.. code-block:: bash

   python -m pip install libhankel


Installing from source
----------------------

You may instead wish to install from source, e.g., to get the very latest version
of the code that is still in development.

1. Download this repository or clone it using:

   .. code-block:: bash

      git clone https://github.com/Letizia97/libhankel.git


2. Open up a terminal (command prompt) and go into the
   ``libhankel`` directory.


3. It is recommended to create a virtual environment by typing for example:

   .. code-block:: bash
      
      python3 -m venv .venv
      . .venv/bin/activate

   If the above fails, you might need to run the following (adjusted to match your Python version):

   .. code-block:: bash

      apt install python3.10-venv

4. Install the `Boost` dependency:

   .. code-block:: bash

      sudo apt update
      sudo apt install libboost-all-dev

5. At this point, we are ready to install LibHankel. 
   Please ensure the virtual environment has been activated through step 3 above, 
   and then type:

   .. code-block:: bash

      python -m pip install .


Please note that this should install all dependencies required (pip dependencies 
are installed automatically). However, in case of issues, the exact dependencies 
are available in the file "requirements_dev.txt".

The docs have multiple Python examples that you can try running at this point.

.. python-installation-end

