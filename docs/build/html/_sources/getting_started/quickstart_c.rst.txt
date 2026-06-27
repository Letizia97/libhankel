Installation for C users
=================================

.. _quickstart-c:


Requirements:

- Meson

- Ninja (recommended)

- A C compiler


Build and install:

.. code-block:: bash

   git clone https://github.com/libhankel
   cd libhankel
   meson setup build
   meson compile -C build
   meson install -C build


After installation, you can compile programs using:

.. code-block:: bash

   pkg-config --cflags --libs libhankel
