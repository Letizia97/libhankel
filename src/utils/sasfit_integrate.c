#include "src/utils/sasfit_integrate.h"

#include <math.h>

#include "../external_libs/utils/tanhsinh.h"
#include "boost_bessel_wrapper.h"
#include "libhankel.h"

typedef double sasfit_func_one_t(double, hankel_inputs *);

typedef struct {
    hankel_inputs *param;
    sasfit_func_one_t *Kernel1D_fct;
} int_cub;

static double Kernel_1D(double x, void *pam) {
    int_cub *cub = (int_cub *)pam;
    return cub->Kernel1D_fct(x, cub->param);
}

double sasfit_integrate_ctm(double int_start, double int_end, sasfit_func_one_t intKern_fct,
                            hankel_inputs *param, double epsrel) {
    int_cub cubstruct;
    double ferr;

    cubstruct.Kernel1D_fct = intKern_fct;
    cubstruct.param = param;

    // nothing to integrate
    if (isfinite(int_start) && isfinite(int_end) && (int_end - int_start) == 0.0) {
        return 0.0;
    }

    return TanhSinhQuad(&Kernel_1D, &cubstruct, int_start, int_end, 7, epsrel, &ferr);
}

// Used in both sasfit_HankelChave and sasfit_qwe
double FrJnu(double r, hankel_inputs *inputs) {
    double Q, nu;
    nu = inputs->other_inputs[0];
    Q = inputs->other_inputs[1];
    return r * bessel_Jnu(nu, Q * r) * inputs->function(r, inputs->f_params);
}
