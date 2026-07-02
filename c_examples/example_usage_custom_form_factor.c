#include <stdio.h>
#include <stddef.h>

#include "libhankel.h"
#include "form_factors.h"

#define ARRAY_LEN 19


int main(void) {

    double nu = 0;

    // Build the context needed for the form factor
    form_factor_ctx ctx_gdab;
    double params_gdab[] = {10.0, 0.5, 1e-4};
    ctx_gdab.params = params_gdab;

    // build the struct required by the chosen strategy
    strategy_params strategy_params_gdab = {
        .n_eval = 250, 
        .eps_rel = 1e-9
    };

    // array where to store the result
    double result[ARRAY_LEN];

    // input of values where to compute the transform
    double r_array_gdab[ARRAY_LEN] = {
        15.,   18.5,  22.,   25.6,  29.1,
        32.7,  36.2,  39.7,  43.3,  46.9,
        50.4,  53.9,  57.5,  61.1,  64.5,
        68.1,  71.6,  75.2,  78.75,             
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

    // print the result
    for (int i = 0; i < ARRAY_LEN; i++) {
        printf("%.15g ", result[i]);
    }
    printf("\n");

    /* Add cleanup */

    return 0;

}

