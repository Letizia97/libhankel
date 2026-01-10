#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_bessel.h>
#include <stdio.h>
#include "include/libhankel.h"

static gsl_integration_workspace *workspace = NULL; 
typedef double  sasfit_func_one_t (double, hankel_inputs *);

typedef struct {
    hankel_inputs *param;
    sasfit_func_one_t *Kernel1D_fct;
} int_cub;

double Kernel_1D(double x, void *pam) {
	hankel_inputs * param;
	int_cub *cub;
	cub = (int_cub *) pam;
	param = (hankel_inputs *) cub->param;
	return (*cub->Kernel1D_fct)(x,param);
}

double sasfit_integrate_ctm(
    double int_start,
    double int_end,
    sasfit_func_one_t intKern_fct,
    hankel_inputs *param,
    int limit,
    double epsabs,
    double epsrel)
{
    workspace = gsl_integration_workspace_alloc(1000);
	double res, errabs; 
	gsl_function F;
    int err;
    int_cub cubstruct;
    double ferr[1];

    cubstruct.Kernel1D_fct=intKern_fct;
    cubstruct.param=param;

    // nothing to integrate
    if ( gsl_finite(int_start) && gsl_finite(int_end) &&
	     (int_end - int_start) == 0.0 )
	{
		return 0.0;
	}

	F.params = param;
	F.function = (double (*) (double, void*)) intKern_fct;

    cubstruct.Kernel1D_fct=intKern_fct;
    cubstruct.param=param;

	gsl_set_error_handler_off();
    err = GSL_SUCCESS;

    res = TanhSinhQuad(&Kernel_1D, &cubstruct, int_start, int_end, 7, epsrel, &ferr[0]);
	return res;
}

double FrJnu(double r, hankel_inputs * inputs) {
    double Q,nu;
    nu = inputs->other_inputs[0];
    Q  = inputs->other_inputs[1];
    return r*gsl_sf_bessel_Jnu(nu,Q*r)*inputs->function(r,inputs->fparams);
}
