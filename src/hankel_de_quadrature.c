
#include "include/libhankel.h"

#include <math.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_bessel.h>

typedef struct
{
	void *fparams; 
	double (* function) (double x, void *);
	double nu;
	double Q;
} params_struct;


double intdeo_FBT(double r, void *FBTparams) {
    /* 
    Function description
    */
    params_struct *FBTparam_struct;
    FBTparam_struct = (params_struct *) FBTparams;
    if (r==0) return 0;
    return r*gsl_sf_bessel_Jnu(FBTparam_struct->nu,FBTparam_struct->Q*r) * FBTparam_struct->function(r,FBTparam_struct->fparams);
}

double hankel_transform_DE_Quadrature(int nu, double (*f)(double, double (*)[50]), double x, double (*fparams)[50], double N_ogata, double eps_nriq){
    /* 
    Computes Hankel transform, using de-quadrature.
    */
    int lenaw=4000;
    int rounded_N;
    double res0, err0, res, err, a, *aw, nv;
    rounded_N = lround(N_ogata);
    nv = fabs(nu);

    params_struct FBTparam_struct;
    FBTparam_struct.fparams=fparams;
    FBTparam_struct.function=f;
    FBTparam_struct.nu=nu;
    FBTparam_struct.Q=x;

    // FIXME: need to find a better variable name for a and aw
    a = gsl_sf_bessel_zero_Jnu(nv, (rounded_N<10?rounded_N:10)) / x;
    aw = (double *)malloc((lenaw)*sizeof(double));

    sasfit_intdeini(lenaw, GSL_DBL_MIN, eps_nriq, aw);
    sasfit_intde(&intdeo_FBT,0, a, aw, &res0, &err0, &FBTparam_struct);

    sasfit_intdeoini(lenaw, GSL_DBL_MIN, eps_nriq, aw);
    sasfit_intdeo(&intdeo_FBT, a, FBTparam_struct.Q, aw, &res, &err, &FBTparam_struct);
    free(aw);

    res += res0;
    err += err0;
    return res;
}