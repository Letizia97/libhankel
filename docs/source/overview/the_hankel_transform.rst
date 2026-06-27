.. _hankel-intro:

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


.. _hankel-transform-strategies:

Hankel Transform Strategies
---------------------------

LibHankel includes various implementations of the Hankel transform, corresponding to 
different ways of approximating the Hankel transform integral.

The following table lists the strategies that are available through this library. 
The column named as "index" in the table essentially refers to the numbering used in SASfit.


.. table:: Hankel transform strategies available in LibHankel.

    +--------------------------+---------------------+---------------------------------------+
    | index (as in SASfit)     | Strategy Name       | Type                                  |
    +==========================+=====================+=======================================+
    | 0                        | DE_Ooura            | Double-exponential quadrature         |
    +--------------------------+---------------------+---------------------------------------+
    | 1                        | DE_Ogata            | Double-exponential quadrature         |
    +--------------------------+---------------------+---------------------------------------+
    | 6 to 11                  | DHT_6 to DHT_11     | Digital filters                       |
    +--------------------------+---------------------+---------------------------------------+
    | 12                       | QWE_Key             | QWE with continued fraction expansion |
    +--------------------------+---------------------+---------------------------------------+
    | 13                       | QWE_Chave           | QWE with Shanks transformation        |
    +--------------------------+---------------------+---------------------------------------+


Strategies 6 to 11 in SASfit, which in this library are referred to as "DHT_6", "DHT_7" ... 
"DHT_11",  are fixed-grid digital filter strategies. They are fast and perform well on simple form
factors but their accuracy decreases on more complex or oscillatory ones. They are not iterative, 
rely on pre-computed constants, and do not allow the user to control the "accuracy" of the solution 
in any way. 

All other strategies instead are iterative, and some of them require the user to supply a tolerance 
parameter (called ``eps_rel``) which is used to determine when to stop interating. Specifically, 
iterations are stopped when the difference in the solution of two consecutive iterations is less than 
the tolerance. For this reasons, strategies that require eps_rel generally perform better than the DHT 
ones on complex form factors, though they are more slow. 

