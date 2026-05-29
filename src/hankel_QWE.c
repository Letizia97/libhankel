#include "libhankel.h"

#include <math.h>
#include <stdio.h>
#include <float.h>

#include "../external_libs/qwe_Key/qwe_Key.h"
#include "../external_libs/qwe_Chave/qwe_Chave.h"

/*
This file contains functions corresponding to strategies 12 and 13 in SASfit. 
They have been grouped together as these are QWE transforms, that is Quadrature With Extrapolation.
Specifically:
- HANKEL_QWE 12
- HANKEL_CHAVE 13
*/


/** 
 * @brief Computes Hankel transform using the Quadrature With Extrapolation method by Key.
 * @note Corresponds to strategy 12 in SASfit.
 * 
 * @param nu         order of bessel function - must be 0 or 1
 * @param f          pointer to form factor function
 * @param x          value at which to compute the transform
 * @param f_params    params for form factor
 * @param output     pointer to var containing output from transform 
 * @param n_eval     integer indicating number of function evaluations (N_ogata in SASfit)
 * @param eps_rel    relative error allowed e.g. 1e-9 (eps_nriq in SASfit)
 */ 
double hankel_transform_QWE_Key(
    int nu, 
    form_factor_f f, 
    const double x,
    double *f_params,
    size_t n_params,
    double * output,
    int n_eval, 
    double eps_rel,
    void *user_data
) {

    int status;
    status = qwe_Key(
        nu, 
        f, 
        x, 
        f_params,
        output,
        lround(n_eval),
        eps_rel*10,
        DBL_MIN);
        
    return status;
}

/** 
 * @brief Computes Hankel transform using the Quadrature With Extrapolation method by Chave.
 * @note Corresponds to strategy 13 in SASfit.
 * 
 * @param nu         order of bessel function - must be 0 or 1
 * @param f          pointer to form factor function
 * @param x          value at which to compute the transform
 * @param f_params    params for form factor
 * @param output     pointer to var containing output from transform 
 * @param n_eval     integer indicating number of function evaluations (N_ogata in SASfit)
 * @param eps_rel    relative error allowed e.g. 1e-9 (eps_nriq in SASfit)
 */ 
double hankel_transform_QWE_Chave(
    int nu, 
    form_factor_f f, 
    const double x,
    double *f_params,
    size_t n_params,
    double * output,
    int n_eval, 
    double eps_rel,
    void *user_data) {

    int status;
    status = qwe_Chave(
        nu, 
        f,
        x, 
        f_params,
        output,
        lround(n_eval),
        eps_rel*10,
        DBL_MIN);

    return status;
}
