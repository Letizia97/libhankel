#include "include/libhankel.h"

#include <math.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_bessel.h>


/*
This file contains functions corresponding to strategies 12 and 13 in SASfit. 
They have been grouped together as these are QWE transforms, that is Quadrature With Extrapolation.
Specifically:
- HANKEL_QWE 12
- HANKEL_CHAVE 13

*/

double hankel_transform_QWE_Key(int nu, double (*f)(double, double (*)[50]), double x, double (*fparams)[50], double N_ogata, double eps_nriq){
    /*
    Func descr
    */
    double res;
    res = sasfit_qwe(
        nu, 
        f, 
        x, 
        fparams,
        lround(N_ogata),
        eps_nriq*10,
        DBL_MIN);
    return res;
    }

double hankel_transform_QWE_Chave(int nu, double (*f)(double, double (*)[50]), double x, double (*fparams)[50], double N_ogata, double eps_nriq){
    /*
    Func descr
    */
    double res;
    res = sasfit_HankelChave(
        nu, 
        f,
        x, 
        fparams,
        lround(N_ogata),
        eps_nriq*10,
        DBL_MIN);
    return res;
    }
