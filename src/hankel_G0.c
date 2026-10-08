#include "libhankel.h"
#include "src/utils/sasfit_integrate.h"
#include <math.h>

static double q_times_f(double q, hankel_inputs *inputs) {
    return q * inputs->function(q, inputs->f_params);
}

double hankel_G0(form_factor_f f, void *f_ctx, double epsrel) {
    hankel_inputs inputs = {0};
    inputs.function = f;
    inputs.f_params = f_ctx;
    return sasfit_integrate_ctm(0.0, INFINITY, q_times_f, &inputs, epsrel);
}
