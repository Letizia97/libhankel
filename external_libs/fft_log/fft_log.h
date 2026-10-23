#pragma once

#ifdef __cplusplus
extern "C" {
#endif

double sasfit_fftlog(double nu, double (*f)(double, void *), double x, void *fparams);

#ifdef __cplusplus
}
#endif
