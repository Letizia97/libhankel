#include "libhankel.h"

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

/** 
 * @brief Computes Hankel transform using FBT.
 * @note Corresponds to strategies 2,3,4 in SASfit or HANKEL_FBT0, etc 
 * @note Requires knowledge of the form factor function, as one of the parameters 
 *       passed (f_max) is a starting guess for maximum in form factor. 
 *       Performance is highly dependent on params.
 *       
 * @param nu         order of bessel function - must be 0 or 1
 * @param f          pointer to form factor function
 * @param x          value at which to compute the transform
 * @param f_params    params for form factor
 * @param n_eval     integer indicating number of function evaluations (N_ogata in SASfit)
 * @param f_max      float indicating starting guess for maximum in form factor (h_ogata in SASfit)
 * @param n_method   FBT method to use among 0,1,2 (modified, unmodified, fixed h Ogata)
 */ 
double hankel_transform_FBT(
    int nu, 
    double (*f)(double, double (*)[50]), 
    double x, 
    double (*f_params)[50], 
    int n_method, 
    double n_eval, 
    double f_max) {

    double res; 
    res = compute_hankel_FBT(nu, f, x, f_params, n_method, n_eval, f_max);
    return res;
}
