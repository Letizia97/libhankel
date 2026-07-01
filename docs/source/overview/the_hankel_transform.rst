.. _hankel-intro:


#######################################################
Intro to LibHankel
#######################################################


The Hankel Transform
---------------------------

The Hankel transform, also known as the Fourier–Bessel transform, is widely used in scattering, specifically for:

- calculating scattering functions;

- correcting for multiple scattering;

- simulating spin echo signals (e.g.,  those obtained on Larmor at the ISIS Neutron and Muon Source);

- qualitatively analysing phase contrast images.

Small-angle scattering (SAS) is an experimental technique using X-rays (SAXS) or neutrons (SANS)
to study the structure of materials on the nanoscale (roughly 1–1000 nm). SAS measures how radiation
is scattered at very small angles, which tells us about large-scale structures inside the sample.


.. _why-libhankel:

Why LibHankel
-------------
`SASfit <https://github.com/SASfit/SASfit>`_ is a tool for analysis of SAS data. LibHankel essentially 
re-implements (some of) the Hankel transform algorithms available in SASfit with the aim of rendering them
easily callable from Python and C, as well as facilitating maintenance.

