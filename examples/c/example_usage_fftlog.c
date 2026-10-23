#include <stddef.h>
#include <stdio.h>

#include "form_factors.h"
#include "libhankel.h"

#define ARRAY_LEN 19

int main(void) {

    int nu = 0;

    /* Build the context for the built-in g_dab form factor */
    form_factor_ctx ctx;
    double params[] = {10.0, 0.5, 1e-4};
    ctx.params = params;

    /* FFTLog takes no strategy parameters; pass a zero-initialised struct */
    strategy_params sp = {0};

    double result[ARRAY_LEN];

    double q_arr[ARRAY_LEN] = {
        15.,  18.5, 22.,  25.6, 29.1, 32.7, 36.2, 39.7, 43.3,  46.9,
        50.4, 53.9, 57.5, 61.1, 64.5, 68.1, 71.6, 75.2, 78.75,
    };

    hankel_transform(nu, form_factor_g_dab, q_arr, ARRAY_LEN, (void *)&ctx, result,
                     "FFTLog", sp);

    for (int i = 0; i < ARRAY_LEN; i++) {
        printf("%.15g ", result[i]);
    }

    return 0;
}
