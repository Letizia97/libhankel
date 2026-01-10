#include "include/libhankel.h"

#include <math.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_bessel.h>

/*
This file contains functions corresponding to strategies 2-4 in SASfit. 
Specifically:
- HANKEL_FBT0 2
- HANKEL_FBT1 3
- HANKEL_FBT2 4
*/

double hankel_transform_FBT(int nu, double (*f)(double, double (*)[50]), double x, double (*fparams)[50], int n_method, double N_ogata, double h_ogata){
    /* 
    Computes Hankel transform, using FBT.
    Requires additional params - specifically number of function evaluations and starting 
    guess for maximum in form factor - and performance is highly dependent on params.
    Requires knowledge of the form factor function to properly set params. 
    Corresponds to strategies 2,3,4 in SASfit or HANKEL_FBT0, etc (specify through n_method 0,1,2 in inputs)
    Receives:
        nu         order of bessel function
        *f         pointer to form factor function
        x          value at which to compute the transform
        *fparams   params for form factor
        N_ogata    integer indicating number of function evaluations
        h_ogata    float indicating starting guess for maximum in form factor 
        n_method   FBT method to use among 0,1,2 (modified, unmodified, fixed h Ogata)
    */ 
    double res; 
    res = compute_hankel_FBT(nu, f, x, fparams, n_method, N_ogata, h_ogata);
    return res;

};
