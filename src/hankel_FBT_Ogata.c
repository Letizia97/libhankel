#include "include/libhankel.h"

#include <math.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_bessel.h>


/*
This file contains functions corresponding to strategies 1-4 in SASfit. 
They have been grouped together as written by Ogata originally.
Specifically:
- HANKEL_OGATA_2005
- HANKEL_FBT0
- HANKEL_FBT1
- HANKEL_FBT2

The first is implemented through function "hankel_transform_DE_Ogata".
The last 3 are implemented through the same function "hankel_transform_FBT".
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


double DEtransform(double t) {
    /* 
    Function descr
    */ 
    return t * tanh(M_PI_2 * sinh(t));
}

double deriv_DEtransform(double t){
    /* 
    Function descr
    */ 
    double res;
    res = 1. / cosh(M_PI_2 * sinh(t));
    res = M_PI_2 * t * cosh(t) * res * res + tanh(M_PI_2 * sinh(t));
    return res;
}

double hankel_transform_DE_Ogata(int nu, double (*f)(double, double (*)[50]), double x, double (*fparams)[50], double N_ogata, double h_ogata){
    /* 
    Function descr
    Corresponds to strategy 1 in SASfit or HANKEL_OGATA_2005.
    */ 
    double res, zeros_PI, phi_dot, y_k, w_nv_k, J_nv, sum, nv;
    int i;
    sum = 0.0;
    nv = abs(nu);

    for (i=1; i<=N_ogata; i++) {
        zeros_PI = gsl_sf_bessel_zero_Jnu(nv,i) / M_PI;
        phi_dot = deriv_DEtransform(h_ogata * zeros_PI);
        y_k = DEtransform(h_ogata * zeros_PI) * M_PI / h_ogata;
        w_nv_k = 2./ (gsl_pow_2(M_PI * gsl_sf_bessel_Jnu(nv + 1, zeros_PI * M_PI)) * zeros_PI);
        J_nv = gsl_sf_bessel_Jnu(nv, y_k);
        res =  w_nv_k * y_k * (*f)(y_k / x, fparams) * J_nv * phi_dot;
        sum = sum + res;
    }

    sum = M_PI / gsl_pow_2(x) * sum;
    if (nv==0) {
        res = sum;
    } else {
        res = sum * pow(GSL_SIGN(nu),nu);
    }
    return res;

};

