#include <stddef.h>
#include <stdio.h>

#include "form_factors.h"
#include "libhankel.h"

int main(void) {

    /* Build the context for the built-in g_dab form factor */
    form_factor_ctx ctx;
    double params[] = {10.0, 0.5, 1e-4};
    ctx.params = params;

    double result = hankel_G0(form_factor_g_dab, (void *)&ctx, 1e-6);

    printf("G(0) = %.15g\n", result);

    return 0;
}
