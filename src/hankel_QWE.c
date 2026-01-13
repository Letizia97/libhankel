#include "include/libhankel.h"

#include <math.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_bessel.h>

#include "external_libs_based_work/qwe/qwe.h"

/*
This file contains functions corresponding to strategies 12 and 13 in SASfit. 
They have been grouped together as these are QWE transforms, that is Quadrature With Extrapolation.
Specifically:
- HANKEL_QWE 12
- HANKEL_CHAVE 13
*/

double hankel_transform_QWE_Key(
    int nu, 
    double (*f)(double, double (*)[50]), 
    double x, 
    double (*fparams)[50], 
    double n_eval, 
    double eps_rel) {
    /* 
    Computes Hankel transform using the Quadrature With Extrapolation method by Key.
    Corresponds to strategy 12 in SASfit.
    Receives:
        nu         order of bessel function (typically 0)
        *f         pointer to form factor function
        x          value at which to compute the transform
        *fparams   params for form factor
        n_eval     integer indicating number of function evaluations (N_ogata in SASfit)
        eps_rel    relative error allowed e.g. 1e-9 (eps_nriq in SASfit)
    */ 
    double res;
    res = sasfit_qwe(
        nu, 
        f, 
        x, 
        fparams,
        lround(n_eval),
        eps_rel*10,
        DBL_MIN);
    return res;
    }

double hankel_transform_QWE_Chave(
    int nu, 
    double (*f)(double, double (*)[50]), 
    double x, 
    double (*fparams)[50], 
    double n_eval, 
    double eps_rel) {
    /* 
    Computes Hankel transform using the Quadrature With Extrapolation method by Key.
    Corresponds to strategy 13 in SASfit.
    Receives:
        nu         order of bessel function (typically 0)
        *f         pointer to form factor function
        x          value at which to compute the transform
        *fparams   params for form factor
        n_eval     integer indicating number of function evaluations (N_ogata in SASfit)
        eps_rel    relative error allowed e.g. 1e-9 (eps_nriq in SASfit)
    */ 
    double res;
    res = sasfit_HankelChave(
        nu, 
        f,
        x, 
        fparams,
        lround(n_eval),
        eps_rel*10,
        DBL_MIN);
    return res;
    }
