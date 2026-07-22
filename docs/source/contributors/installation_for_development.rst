
.. _installation_for_development:

Installation for development
==============================

Installation for a contributor is very similar to that for a C user.
As mentioned in other parts of the docs, LibHankel requires:

- Meson >= 1.4.0
- Ninja
- A C compiler (e.g. gcc)

If installing Meson with ``pip``, it is recommended to use a Python virtual environment to avoid
modifying system Python packages.

.. code-block:: bash

   python3 -m venv .venv
   source .venv/bin/activate
   python -m pip install meson==1.4.0


LibHankel also requires Boost development package as a dependency. This can be installed with:

.. code-block:: bash

   sudo apt update
   sudo apt install libboost-all-dev

To build and install LibHankel, together with the test suite, please use:

.. code-block:: bash

   git clone --recurse-submodules https://github.com/Letizia97/libhankel.git
   cd libhankel
   meson setup build -Dtests=true
   meson compile -C build
   meson install -C build

Please note that using ``git clone`` and ``meson setup build`` without the suggested flags will not 
include the tests in the installation.
