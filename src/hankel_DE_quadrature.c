#include "libhankel.h"

#include <math.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_bessel.h>

#include "external_libs_based_work/DE-quadrature/intde.h"

/*
This file contains only 2 Hankel strategies, both DE quadrature algorithms, corresponding to:
    - HANKEL_OOURA_DEO , i.e. strategy 0 in SASfit
    - OGATA_2005 , i.e. strategy 1 in SASfit
*/


/**
 * @brief Struct of parameters to be used in DE hankel functions.
 */
typedef struct {
	void *fparams;  /**< parameters for the supplied function */
	double (* function) (double, double (*)[50]);  /**< function to integrate */
	double nu;  /**< order of the Bessel function */
	double Q;   /**< radial Fourier variable, i.e. conj wavenumber to radius r */
} params_struct;


/** 
 * @brief Auxiliary function that computes the value of the 
 *        Hankel‑transform integrand at the current radius r.
 *        Computes a product between the radius r, the Bessel 
 *        function of the first kind of order nu, and a function
 *        supplied through FBTparams.
 * 
 * @param r          radius
 * @param FBTparams  pointer to a struct containing nu, Q and function to integrate
 */
double intdeo_FBT(double r, void *FBTparams) {
    params_struct *FBTparam_struct;
    FBTparam_struct = (params_struct *) FBTparams;
    if (r==0) return 0;

    double nu = FBTparam_struct->nu;
    double Q  = FBTparam_struct->Q;

    double bessel = gsl_sf_bessel_Jnu(nu, Q * r);
    double fval   = FBTparam_struct->function(r, FBTparam_struct->fparams);

    return r * bessel * fval;
}

/** 
 * @brief Auxiliary function computing the 
 *        double‑exponential (DE) / tanh–sinh transform.
 */
double DEtransform(double t) {
    return t * tanh(M_PI_2 * sinh(t));
}

/** 
 * @brief Auxiliary function that computes the derivative of a 
 *        “double‑exponential (DE) change of variables”, that is 
 *        commonly used in numerical integration (quadrature).
 *        Specifically, this is the tanh–sinh (double exponential) 
 *        transformation.
 */
double deriv_DEtransform(double t){
    double res;
    res = 1. / cosh(M_PI_2 * sinh(t));
    res = M_PI_2 * t * cosh(t) * res * res + tanh(M_PI_2 * sinh(t));
    return res;
}

/** 
 * @brief Computes Hankel transform, using de-quadrature.
 * @note Corresponds to strategy 0 in SASfit, or HANKEL_OOURA_DEO.
 *        
 * @param nu         order of bessel function - must be 0 or 1
 * @param f          pointer to form factor function
 * @param x          value at which to compute the transform
 * @param fparams    params for form factor
 * @param n_eval     integer indicating number of function evaluations (N_ogata in SASfit)
 * @param eps_rel    relative error allowed e.g. 1e-9 (eps_nriq in SASfit)
 */ 
double hankel_transform_DE_Quadrature(
    int nu, 
    double (*f)(double, double (*)[50]), 
    double x, 
    double (*fparams)[50], 
    double n_eval, 
    double eps_rel) {

    int lenaw = 4000;
    int rounded_N;
    double res0, err0, res, err, a, *aw, nv;
    rounded_N = lround(n_eval);
    nv = abs(nu);

    params_struct FBTparam_struct;
    FBTparam_struct.fparams=fparams;
    FBTparam_struct.function=f;
    FBTparam_struct.nu=nu;
    FBTparam_struct.Q=x;

    // FIXME: need to find a better variable name for a and aw
    a = gsl_sf_bessel_zero_Jnu(nv, (rounded_N<10?rounded_N:10)) / x;
    aw = (double *)malloc((lenaw)*sizeof(double));

    sasfit_intdeini(lenaw, GSL_DBL_MIN, eps_rel, aw);
    sasfit_intde(&intdeo_FBT,0, a, aw, &res0, &err0, &FBTparam_struct);

    sasfit_intdeoini(lenaw, GSL_DBL_MIN, eps_rel, aw);
    sasfit_intdeo(&intdeo_FBT, a, FBTparam_struct.Q, aw, &res, &err, &FBTparam_struct);
    free(aw);

    res += res0;
    err += err0;
    return res;
}

/** 
 * @brief Computes Hankel transform, using de-quadrature.
 * @note Corresponds to strategy 1 in SASfit or HANKEL_OGATA_2005.
 * 
 * @param nu         order of bessel function - must be 0 or 1
 * @param f          pointer to form factor function
 * @param x          value at which to compute the transform
 * @param fparams    params for form factor
 * @param n_eval     integer indicating number of function evaluations (N_ogata in SASfit)
 * @param f_max      float indicating starting guess for max in form factor (h_ogata in SASfit)
 */ 
double hankel_transform_DE_Ogata(
    int nu, 
    double (*f)(double, double (*)[50]), 
    double x, 
    double (*fparams)[50], 
    double n_eval, 
    double f_max) {

    double res, zeros_PI, phi_dot, y_k, w_nv_k, J_nv, sum, nv;
    int i;
    sum = 0.0;
    nv = abs(nu);

    for (i=1; i<=n_eval; i++) {
        zeros_PI = gsl_sf_bessel_zero_Jnu(nv,i) / M_PI;
        phi_dot = deriv_DEtransform(f_max * zeros_PI);
        y_k = DEtransform(f_max * zeros_PI) * M_PI / f_max;
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

}

