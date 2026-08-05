
#include "libhankel.h"

// Standard library headers
#include <math.h>
#include <stdio.h>

// Project / local headers
#include "../src/utils/strateg10_const.h"
#include "../src/utils/strateg11_const.h"
#include "../src/utils/strateg6_const.h"
#include "../src/utils/strateg7_const.h"
#include "../src/utils/strateg8_const.h"
#include "../src/utils/strateg9_const.h"
#include "../src/utils/validate_x.h"

/*
This file contains functions corresponding to strategies 6-11 in SASfit
(see inline comments for specific names).

They have been grouped together under one function, as they are very similar,
and changing the n_strategy parameter allows to switch between them
*/

/*
Strategies 6 and 7 place their nodes at 10^(a + i*s), which depends only on the
filter, not on x or on the form factor. Evaluating those powers on every call
costs more than the rest of the loop put together, so each filter fills its
table on first use and indexes it from then on.
*/
#define GS_MAX_NODES 140

typedef struct {
    double node[GS_MAX_NODES];
    int filled;
} gs_nodes;

static const double *gs_node_table(gs_nodes *table, double a, double s, unsigned int n) {
    unsigned int i;

    if (!table->filled) {
        for (i = 0; i < n; i++) {
            table->node[i] = pow(10.0E0, (a + i * s));
        }
        table->filled = 1;
    }
    return table->node;
}

int hankel_transform_DHT(int nu, form_factor_f f, const double x, void *f_ctx, double *output,
                         int n_strategy) {

    static gs_nodes gs_J0, gs_J1, gs_J0_fast, gs_J1_fast;

    double res = 0;
    double lambda;
    const double *node;
    unsigned int i;
    unsigned int ind;

    if (!(nu == 0 || nu == 1)) {
        fprintf(stderr, "nu needs to be 0 or 1 in order to use the selected strategy\n");
        return -1;
    }

    ind = nu + 1;

    if (!(n_strategy >= 6 && n_strategy <= 11)) {
        fprintf(stderr, "Strategy number must be integer between 6 and 11\n");
        return -2;
    }

    /* Every filter below samples the form factor at lambda = node / x and
     * weights it by 1 / x, so a non-positive or non-finite x would return
     * NaN - or, for x < 0, a plausible-looking finite number - as a success. */
    int status = validate_x(x);
    if (status != 0) {
        return status;
    }

    switch (n_strategy) {

    case 6: {
        // HANKEL_GUPTASARMA_97
        if (nu == 0) {
            node = gs_node_table(&gs_J0, aJ0, sJ0, 120);
            for (i = 0; i < 120; i++) {
                lambda = node[i] / x;
                res = res + (*f)(lambda, f_ctx) * lambda * WJ0[i] / x;
            }
        } else {
            node = gs_node_table(&gs_J1, aJ1, sJ1, 140);
            for (i = 0; i < 140; i++) {
                lambda = node[i] / x;
                res = res + (*f)(lambda, f_ctx) * lambda * WJ1[i] / x;
            }
        }
        break;
    }
    case 7: {
        // HANKEL_GUPTASARMA_97_FAST
        if (nu == 0) {
            node = gs_node_table(&gs_J0_fast, aJ0Fast, sJ0Fast, 61);
            for (i = 0; i < 61; i++) {
                lambda = node[i] / x;
                res = res + (*f)(lambda, f_ctx) * lambda * WJ0Fast[i] / x;
            }
        } else {
            node = gs_node_table(&gs_J1_fast, aJ1Fast, sJ1Fast, 47);
            for (i = 0; i < 47; i++) {
                lambda = node[i] / x;
                res = res + (*f)(lambda, f_ctx) * lambda * WJ1Fast[i] / x;
            }
        }
        break;
    }
    case 8: {
        // HANKEL_KEY_51
        for (i = 0; i < 51; i++) {
            lambda = KK51Hankel[i][0] / x;
            res = res + (*f)(lambda, f_ctx) * lambda * KK51Hankel[i][ind] / x;
        }
        break;
    }
    case 9: {
        // HANKEL_KEY_101
        for (i = 0; i < 101; i++) {
            lambda = KK101Hankel[i][0] / x;
            res = res + (*f)(lambda, f_ctx) * lambda * KK101Hankel[i][ind] / x;
        }
        break;
    }
    case 10: {
        // HANKEL_KEY_201
        for (i = 0; i < 201; i++) {
            lambda = KK201Hankel[i][0] / x;
            res = res + (*f)(lambda, f_ctx) * lambda * KK201Hankel[i][ind] / x;
        }
        break;
    }
    case 11: {
        // HANKEL_ANDERSON_801
        for (i = 0; i < 801; i++) {
            lambda = WA801Hankel[i][0] / x;
            res = res + (*f)(lambda, f_ctx) * lambda * WA801Hankel[i][ind] / x;
        }
        break;
    }
    default:
        /* Unreachable: validated above */
        return -2;
    }
    *output = res;
    return 0;
}
