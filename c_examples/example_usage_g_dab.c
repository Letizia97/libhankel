#include <stdio.h>
#include <stddef.h>

#include "libhankel.h"
#include "form_factors.h"

#define ARRAY_LEN 19


int main(void) {

    double nu = 0;
    form_factor_ctx ctx_gdab;
    double params_gdab[] = {10.0, 0.5, 1e-4};
    ctx_gdab.params = params_gdab;

    strategy_params strategy_params_gdab = {
        .n_eval = 250, 
        .eps_rel = 1e-9
    };

    double result[ARRAY_LEN];

    double r_array_gdab[ARRAY_LEN] = {
        15.,          18.54166667,  22.08333333,  25.625,       29.16666667,
        32.70833333,  36.25,        39.79166667,  43.33333333,  46.875,
        50.41666667,  53.95833333,  57.5,         61.04166667,  64.58333333,
        68.125,       71.66666667,  75.20833333,  78.75,             
    };

    hankel_transform(
        nu, 
        form_factor_g_dab, 
        r_array_gdab, 
        ARRAY_LEN,
        (void *)&ctx_gdab,
        result, 
        "QWE_Key", 
        strategy_params_gdab
    );

    for (int i = 0; i < ARRAY_LEN; i++) {
        printf("%.15g ", result[i]);
    }
    printf("\n");

    /* Add cleanup */

    return 0;

}


// To run (will need to add these to the docs)
// meson compile -C build
// meson install -C build
// sudo ldconfig (refreshes the system library cache)
// gcc  c_examples/example_usage_g_dab.c -llibhankel -o example_usage_g_dab
// ./example_usage_g_dab
