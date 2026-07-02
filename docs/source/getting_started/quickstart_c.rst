.. _quickstart-c:

Installation for C users
=================================

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


After installation, you can compile programs using:

.. code-block:: bash

   pkg-config --cflags --libs libhankel
