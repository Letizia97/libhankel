#include "libhankel.h"

#include <math.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_bessel.h>

#include "../external_libs/DE-quadrature/intde.h"

/*
This file contains only 2 Hankel strategies, both DE quadrature algorithms, corresponding to:
    - HANKEL_OOURA_DEO , i.e. strategy 0 in SASfit
    - OGATA_2005 , i.e. strategy 1 in SASfit
*/


/**
 * @brief Struct of parameters to be used in DE hankel functions.
 */
typedef struct {
	void *f_params;  /**< parameters for the supplied function */
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
    double fval   = FBTparam_struct->function(r, FBTparam_struct->f_params);

    return r * bessel * fval;
}

/** 
 * @brief Auxiliary function computing the 
 *        double‑exponential (DE) / tanh–sinh transform.
 */
double DEtransform(double t) {
    double s = sinh(t);
    double a = M_PI_2 * s;
    return t * tanh(a);
}

/** 
 * @brief Auxiliary function that computes the derivative of a 
 *        “double‑exponential (DE) change of variables”, that is 
 *        commonly used in numerical integration (quadrature).
 *        Specifically, this is the tanh–sinh (double exponential) 
 *        transformation.
 */
double deriv_DEtransform(double t){
    double sh = sinh(t);
    double ch = cosh(t);
    double A  = M_PI_2 * sh;
    double secH = 1.0 / cosh(A);  // sech(A)

    return M_PI_2 * t * ch * (secH * secH) + tanh(A);
}

/** 
 * @brief Computes Hankel transform, using de-quadrature.
 * @note Corresponds to strategy 0 in SASfit, or HANKEL_OOURA_DEO.
 *        
 * @param nu         order of bessel function - must be 0 or 1
 * @param f          pointer to form factor function
 * @param x          value at which to compute the transform
 * @param f_params    params for form factor
 * @param output     pointer to var containing output from transform 
 * @param n_eval     integer indicating number of function evaluations (N_ogata in SASfit)
 * @param eps_rel    relative error allowed e.g. 1e-9 (eps_nriq in SASfit)
 */ 
double hankel_transform_DE_Ooura(
    int nu, 
    double (*f)(double, double (*)[50]), 
    double x, 
    double (*f_params)[50], 
    double *output,
    int n_eval, 
    double eps_rel) {

    int workspace_len = 4000;
    int rounded_n_eval, status;
    double res0, err0, res, err;
    rounded_n_eval = lround(n_eval);

    params_struct FBTparam_struct;
    FBTparam_struct.f_params = f_params;
    FBTparam_struct.function = f;
    FBTparam_struct.nu = nu;
    FBTparam_struct.Q = x;

    // chooses which zero index to request, caps it at 10
    double zero_index = rounded_n_eval < 10 ? rounded_n_eval : 10;

    // compute zero through bessel function and scale it
    double scaled_zero = gsl_sf_bessel_zero_Jnu(nu, zero_index) / x;

    // allocates an array of doubles and returns a pointer to it
    double *workspace = malloc(workspace_len * sizeof *workspace);

    // precompute nodes & weights for DE integration on a finite interval [a, b]
    sasfit_intdeini(
        workspace_len, 
        GSL_DBL_MIN, 
        eps_rel, 
        workspace
    );
    // compute integral using DE quadrature with weights created by sasfit_intdeini
    sasfit_intde(
        &intdeo_FBT, 
        0, 
        scaled_zero, 
        workspace, 
        &res0, 
        &err0, 
        &FBTparam_struct
    );
    // precompute nodes/weights for oscillatory integrals
    // e.g. f(x) cos(omega x) , over [a, inf]
    sasfit_intdeoini(
        workspace_len, 
        GSL_DBL_MIN, 
        eps_rel, 
        workspace
    );
    // evaluate oscillatory integrals using the table built by intdeoini.
    sasfit_intdeo(
        &intdeo_FBT, 
        scaled_zero, 
        FBTparam_struct.Q, 
        workspace, 
        &res, 
        &err, 
        &FBTparam_struct
    );

    free(workspace);
    res += res0;

    // estimate of the numerical integration error for the computed integral 
    // FIXME: need to figure out what to do with this
    err += err0; 

    *output = res;
	return 0;
}

/** 
 * @brief Computes Hankel transform, using de-quadrature.
 * @note Corresponds to strategy 1 in SASfit or HANKEL_OGATA_2005.
 * 
 * @param nu         order of bessel function - must be 0 or 1
 * @param f          pointer to form factor function
 * @param x          value at which to compute the transform
 * @param f_params    params for form factor
 * @param output     pointer to var containing output from transform 
 * @param n_eval     integer indicating number of function evaluations (N_ogata in SASfit)
 * @param f_max      float indicating starting guess for max in form factor (h_ogata in SASfit)
 */ 
double hankel_transform_DE_Ogata(
    int nu, 
    double (*f)(double, double (*)[50]), 
    double x, 
    double (*f_params)[50], 
    double *output,
    int n_eval, 
    double f_max) {

    double sum;
    int status;
    sum = 0.0;

    for (int i = 1; i <= n_eval; i++) {

        /* ---- Get Bessel zero α_{ν,i} scaled by π ---- */
        double zero_i      = gsl_sf_bessel_zero_Jnu(nu, i);
        double zero_scaled = zero_i / M_PI;

        /* ---- Apply DE transform & its derivative ---- */
        double t           = f_max * zero_scaled;
        double phi         = DEtransform(t);
        double phi_prime   = deriv_DEtransform(t);   /* Jacobian */

        /* ---- Map DE node into actual integration node y_k ---- */
        double y_k         = phi * (M_PI / f_max);

        /* ---- Precompute Bessel factors ---- */
        double Jnu_yk      = gsl_sf_bessel_Jnu(nu, y_k);
        double Jnu1_zero   = gsl_sf_bessel_Jnu(nu + 1, zero_i);

        /* ---- Quadrature weight for α_{ν,i} ---- */
        double denom       = M_PI * Jnu1_zero;
        double weight      = 2.0 / ( (denom * denom) * zero_scaled );

        /* ---- Evaluate integrand at scaled location ---- */
        double f_val       = (*f)(y_k / x, f_params);

        /* ---- Assemble quadrature contribution ---- */
        double term        = weight * y_k * f_val * Jnu_yk * phi_prime;

        sum += term;
    }

    // Apply normalization from the change of variables t = x * r
    double scaled_sum = (M_PI / (x * x)) * sum;

    // Parity correction for integer Bessel order (H_-ν = (-1)^ν H_ν)
    double parity = (nu == 0) ? 1.0 : pow((double)GSL_SIGN(nu), nu);

    // Final result
    double res = scaled_sum * parity;

    *output = res;
	return 0;
}

