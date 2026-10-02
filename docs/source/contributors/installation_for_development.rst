
.. _installation_for_development:

Installation for development
==============================

Installation for a contributor is very similar to that for a C user.
As mentioned in other parts of the docs, LibHankel requires:

- Meson >= 1.4.0
- Ninja
- A C compiler and a C++ compiler (e.g. gcc/g++ on Linux, Apple clang on macOS, MinGW-w64 GCC on Windows)

Due to the required Meson version being no less than 1.4.0, it will be necessary
to install Meson through ``pip`` (installing through ``apt update`` most likely won't work).
It is recommended to use a Python virtual environment to avoid modifying system Python
packages:

.. code-block:: bash

   python3 -m venv .venv
   source .venv/bin/activate
   python -m pip install meson==1.4.0


Ninja is generally installed together with Meson. If it isn't, it can be installed with:

.. code-block:: bash

   sudo apt install ninja-build      # Debian / Ubuntu
   brew install ninja                # macOS
   pip install ninja                 # Windows


LibHankel also requires Boost development headers as a dependency:

.. code-block:: bash

   sudo apt update                   # Debian / Ubuntu
   sudo apt install libboost-dev

   brew install boost                # macOS

   vcpkg install boost-math          # Windows (vcpkg)

To build and install LibHankel, together with the test suite, please use:

.. code-block:: bash

   git clone --recurse-submodules https://github.com/Letizia97/libhankel.git
   cd libhankel
   meson setup build -Dtests=true
   meson compile -C build
   meson install -C build

Please note that using ``git clone`` and ``meson setup build`` without the suggested flags will not
include the tests in the installation.

.. note::

   On Windows, use `MSYS2 <https://www.msys2.org/>`_ to install MinGW-w64 and make sure its
   ``bin`` directory is on your ``PATH`` before running Meson.  The vcpkg Boost headers are found
   automatically; if you install them elsewhere, pass ``CPPFLAGS=-I/path/to/boost`` to Meson.
