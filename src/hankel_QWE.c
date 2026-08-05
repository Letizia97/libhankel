#include "../external_libs/qwe_Chave/qwe_Chave.h"
#include "../external_libs/qwe_Key/qwe_Key.h"
#include "libhankel.h"

// Standard library headers
#include <float.h>
#include <math.h>

// Project / local headers
#include "../src/utils/validate_x.h"

/*
This file contains functions corresponding to strategies 12 and 13 in SASfit.
They have been grouped together as these are QWE transforms, that is Quadrature
With Extrapolation. Specifically:
- HANKEL_QWE 12
- HANKEL_CHAVE 13
*/

/* Both routines below place their quadrature intervals at the Bessel zeros
 * divided by x, so a non-positive or non-finite x would return DBL_MAX, Inf or
 * - for x < 0 - a plausible-looking finite number, all reported as a success. */

int hankel_transform_QWE_Key(int nu, form_factor_f f, const double x, void *f_ctx, double *output,
                             int n_eval, double eps_rel) {

    int status;

    status = validate_x(x);
    if (status != 0) {
        return status;
    }

    status = qwe_Key(nu, f, x, f_ctx, output, lround(n_eval), eps_rel * 10, DBL_MIN);

    return status;
}

int hankel_transform_QWE_Chave(int nu, form_factor_f f, const double x, void *f_ctx, double *output,
                               int n_eval, double eps_rel) {
    int status;

    status = validate_x(x);
    if (status != 0) {
        return status;
    }

    status = qwe_Chave(nu, f, x, f_ctx, output, lround(n_eval), eps_rel * 10, DBL_MIN);

    return status;
}
