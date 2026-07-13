#include <stddef.h>
#include <stdio.h>

#include "form_factors.h"
#include "libhankel.h"
#include <math.h>

#define ARRAY_LEN 19

#include "../src/utils/pow_functions.h"
#include "../src/utils/sf_functions.h"

/**
 * Computes the g_dab form factor.
 *
 * Parameters:
 *  q          value at which to compute the transform.
 *  f_ctx      pointer to struct containing inputs for form factor.
 *             f_ctx must contain a pointer to an array with the
 *             following parameters, in this order:
 *              -   xi: correlation length (e.g., 10.0)
 *              -   H : Hurst exponent (e.g., 0.5)
 *              -   ETA : scattering length density contrast (e.g., 1e-4)
 */
double g_dab(double q, void *f_ctx) {

    // Cast the void pointer `f_ctx` to its real type (`form_factor_ctx *`)
    // This is required because void* cannot be dereferenced directly
    form_factor_ctx *ctx = (form_factor_ctx *)f_ctx;
    double *params = ctx->params;

    double XI = params[0];
    double H = params[1];
    double ETA = params[2];
    double numer, denom;

    double factor1 = pow3(2 * XI);
    double factor2 = sf_poch(H, 1.5);

    numer = pow2(factor1 * factor2 * ETA) * pow3(M_PI);
    denom = pow(1 + pow2(q * XI), 1.5 + H);
    return numer / denom;
}

int main(void) {

    double nu = 0;

    // Build the context needed for the form factor
    form_factor_ctx ctx_gdab;
    double params_gdab[] = {10.0, 0.5, 1e-4};
    ctx_gdab.params = params_gdab;

    // build the struct required by the chosen strategy
    strategy_params strategy_params_gdab = {.n_eval = 250, .eps_rel = 1e-9};

    // array where to store the result
    double result[ARRAY_LEN];

    // input of values where to compute the transform
    double r_array_gdab[ARRAY_LEN] = {
        15.,  18.5, 22.,  25.6, 29.1, 32.7, 36.2, 39.7, 43.3,  46.9,
        50.4, 53.9, 57.5, 61.1, 64.5, 68.1, 71.6, 75.2, 78.75,
    };

    hankel_transform(nu, g_dab, r_array_gdab, ARRAY_LEN, (void *)&ctx_gdab, result, "QWE_Key",
                     strategy_params_gdab);

    // print the result
    for (int i = 0; i < ARRAY_LEN; i++) {
        printf("%.15g ", result[i]);
    }

    return 0;
}
