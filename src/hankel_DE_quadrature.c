#include "include/libhankel.h"

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


typedef struct
{
	void *fparams; 
	double (* function) (double, double (*)[50]); //(double x, void *);
	double nu;
	double Q;
} params_struct;


double intdeo_FBT(double r, void *FBTparams) {
    /* 
    Auxiliary function
    */
    params_struct *FBTparam_struct;
    FBTparam_struct = (params_struct *) FBTparams;
    if (r==0) return 0;
    return r*gsl_sf_bessel_Jnu(FBTparam_struct->nu,FBTparam_struct->Q*r) * FBTparam_struct->function(r,FBTparam_struct->fparams);
}

double DEtransform(double t) {
    /* 
    Auxiliary function
    */ 
    return t * tanh(M_PI_2 * sinh(t));
}

double deriv_DEtransform(double t){
    /* 
    Auxiliary function
    */ 
    double res;
    res = 1. / cosh(M_PI_2 * sinh(t));
    res = M_PI_2 * t * cosh(t) * res * res + tanh(M_PI_2 * sinh(t));
    return res;
}


double hankel_transform_DE_Quadrature(
    int nu, 
    double (*f)(double, double (*)[50]), 
    double x, 
    double (*fparams)[50], 
    double n_eval, 
    double eps_rel) {
    /* 
    Computes Hankel transform, using de-quadrature.
    Corresponds to strategy 0 in SASfit, or HANKEL_OOURA_DEO.

    Receives:
        nu         order of bessel function (typically 0)
        *f         pointer to form factor function
        x          value at which to compute the transform
        *fparams   params for form factor
        n_eval     integer indicating number of function evaluations (N_ogata in SASfit)
        eps_rel    relative error allowed e.g. 1e-9 (eps_nriq in SASfit)
    */ 
    int lenaw=4000;
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


double hankel_transform_DE_Ogata(
    int nu, 
    double (*f)(double, double (*)[50]), 
    double x, 
    double (*fparams)[50], 
    double n_eval, 
    double f_max) {
    /* 
    Function descr
    Corresponds to strategy 1 in SASfit or HANKEL_OGATA_2005.

    Receives:
        nu         order of bessel function (typically 0)
        *f         pointer to form factor function
        x          value at which to compute the transform
        *fparams   params for form factor
        n_eval     integer indicating number of function evaluations (N_ogata in SASfit)
        f_max      float indicating starting guess for max in form factor (h_ogata in SASfit)
    */ 
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

};

