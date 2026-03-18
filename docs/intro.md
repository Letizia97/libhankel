# Introduction {#intro}

The Hankel transform, also known as the Fourier–Bessel transform, is widely used in scattering. The
transform is needed in many areas, e.g. to calculate scattering functions, to correct for multiple scattering,
to simulate spin echo signals like those obtained on Larmor at the ISIS Neutron and Muon Source, as
well as to qualitatively analyse phase contrast images (PCI) taken at the Diamond synchrotron source.

Small-angle scattering (SAS) is an experimental technique using X-rays (SAXS) or neutrons (SANS)
to study the structure of materials on the nanoscale (roughly 1–1000 nm). SAS measures how radiation
is scattered at very small angles, which tells us about large-scale structures inside the sample.

In this library, there are currently ... different implementations of the Hankel transform, correspnding to different ways of approximating the Hankel transform integral.


| index     | Strategy Name       | Type                                  |
|:----------|:-------------------:|---------------------------------------|
| 0         | DE                  | Double-exponential quadrature         |
| 1         | DE                  | Double-exponential quadrature         |
| 12        | QWE_Key             | QWE with continued fraction expansion |
| 13        | QWE_Chave           | QWE with Shanks transformation        |
| ...       | ...                 | ...                                   |
