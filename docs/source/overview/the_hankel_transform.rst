.. _hankel-intro:


#######################################################
Intro to LibHankel
#######################################################


The Hankel Transform
---------------------------

The Hankel transform of order :math:`\nu`, also known as the Fourier–Bessel transform, is defined as

.. math::

    G(x) = \int_0^\infty q \, f(q) \, J_\nu(x q) \, dq

where :math:`f(q)` is the form factor, :math:`J_\nu` is the Bessel function of the first kind of order
:math:`\nu`, and :math:`x` is the transform variable. LibHankel supports :math:`\nu = 0` and :math:`\nu = 1`.

Small-angle scattering (SAS) is an experimental technique using X-rays (SAXS) or neutrons (SANS)
to study the structure of materials on the nanoscale (roughly 1–1000 nm). It measures how radiation
is scattered at very small angles, revealing large-scale structures inside the sample.
The Hankel transform arises naturally in this context because the scattering signal and the
real-space structure are related by a transform of this form. It appears specifically when:

- computing scattering functions from a structural model;
- correcting for multiple scattering;
- simulating spin echo signals (e.g., in SESANS experiments such as those carried out on Larmor at the ISIS Neutron and Muon Source);
- analysing phase contrast images.


.. _why-libhankel:

Why LibHankel
-------------
`SASfit <https://github.com/SASfit/SASfit>`_ is a widely used tool for the analysis of SAS data.
LibHankel re-implements a selection of the Hankel transform algorithms available in SASfit,
with two goals: making them directly callable from C and Python without depending on the full
SASfit installation, and providing a clean, maintainable codebase that is straightforward to
extend with new strategies.

