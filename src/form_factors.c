#include "form_factors.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "../src/utils/pow_functions.h"
#include "../src/utils/sf_functions.h"
#include <stddef.h>
#include <string.h>

double form_factor_g_dab(double q, void *f_ctx) {

    // Cast the void pointer `f_ctx` to its real type (`form_factor_ctx *`)
    // This is required because void* cannot be dereferenced directly
    form_factor_ctx *ctx = (form_factor_ctx *)f_ctx;
    double *params = ctx->params;

    double XI = params[0];
    double H = params[1];
    double ETA = params[2];
    double numer, denom;

    double factor1 = pow3(2*XI);
    double factor2 = sf_poch(H,1.5);

    numer = pow2(factor1 * factor2 * ETA) * pow3(M_PI);
    denom = pow(1 + pow2(q * XI), 1.5 + H);
    return numer / denom;
}


double form_factor_sphere(double q, void *f_ctx) {

    double *params = ((form_factor_ctx *)f_ctx)->params;
    double R = params[0];
    double ETA = params[1];
    double interm, factor;

	if (q * R < 1e-4) {
        factor = 1 - pow2(q*R)/10. + pow4(q*R)/280. - pow6(q*R)/15120.;
		interm = ETA * 4.0/3.0 * M_PI * R * R * R * factor;
	} else {
        factor = sin(q*R) - q*R*cos(q*R);
		interm = ETA * 4.0 * M_PI * factor / pow3(q);
	}
    return pow2(interm);
}


double form_factor_broad_peak(double q, void *f_ctx) {

    double *params = ((form_factor_ctx *)f_ctx)->params;
    double I0 = params[0];
    double XI = params[1];
    double Q0 = params[2];
    double M = params[3];
    double P = params[4];

    double interm = 1.0 + pow(fabs(q-Q0) * XI, M);
	return I0 / pow(interm, P);
}

typedef struct {
    const char *name;
    form_factor_f func;
} entry;

static entry table[] = {
    {"gdab", form_factor_g_dab},
    {"broad_peak", form_factor_broad_peak},
    {"sphere", form_factor_sphere},
    {NULL, NULL}
};

form_factor_f get_form_factor_by_name(const char *name) {
    for (int i = 0; table[i].name != NULL; i++) {
        if (strcmp(name, table[i].name) == 0)
            return table[i].func;
    }
    return NULL;
}
