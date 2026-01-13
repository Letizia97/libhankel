#include "include/libhankel.h"

#include <math.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_bessel.h>

// #include "src/utils/sasfit_integrate.h"

/*
This file contains functions corresponding to strategies 2-4 in SASfit. 
Specifically:
- HANKEL_FBT0 2
- HANKEL_FBT1 3
- HANKEL_FBT2 4
*/

double hankel_transform_FBT(
    int nu, 
    double (*f)(double, double (*)[50]), 
    double x, 
    double (*fparams)[50], 
    int n_method, 
    double n_eval, 
    double f_max) {
    /* 
    Computes Hankel transform using FBT.
    Requires knowledge of the form factor function, as one of the parameters passed (f_max)
    is a starting guess for maximum in form factor.
    Also requires number of function evaluations (n_eval).
    Performance is highly dependent on params.
    Corresponds to strategies 2,3,4 in SASfit or HANKEL_FBT0, etc (specify through n_method 0,1,2 in inputs)
    Receives:
        nu         order of bessel function (typically 0)
        *f         pointer to form factor function
        x          value at which to compute the transform
        *fparams   params for form factor
        n_eval     integer indicating number of function evaluations (N_ogata in SASfit)
        f_max      float indicating starting guess for maximum in form factor (h_ogata in SASfit)
        n_method   FBT method to use among 0,1,2 (modified, unmodified, fixed h Ogata)
    */ 
    double res; 
    res = compute_hankel_FBT(nu, f, x, fparams, n_method, n_eval, f_max);
    return res;

};
