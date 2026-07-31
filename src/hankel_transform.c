#include "libhankel.h"
#include <stddef.h>
#include <stdio.h>
#include <string.h>

int validate_n_eval(strategy_params strategy_params) {
    if (strategy_params.n_eval == 0) {
        fprintf(stderr, "Error: n_eval must be provided and cannot be zero\n");
        return -8;
    }
    return 0;
}

int validate_eps_rel(strategy_params strategy_params) {
    if (strategy_params.eps_rel == 0) {
        fprintf(stderr, "Error: eps_rel must be provided and cannot be zero\n");
        return -9;
    }
    return 0;
}

int validate_f_max(strategy_params strategy_params) {
    if (strategy_params.f_max == 0) {
        fprintf(stderr, "Error: f_max must be provided and cannot be zero\n");
        return -10;
    }
    return 0;
}

/* ------------------------------------------------------------------------- *
 * Strategy names
 *
 * Names are FAMILY[_Author][_ntaps]:
 *
 *   DHT_*   digital linear filters - always fixed-tap, named after the published
 *           filter they implement ("DHT_Key_101" is Key's 101-point filter, not a
 *           refinement of "DHT_Key_51"; the three authors are independent designs
 *           and accuracy is NOT monotonic in the tap count).  The Guptasarma pair
 *           carries no tap count because theirs depends on the order: 120 taps for
 *           J0 but 140 for J1 (61 / 47 for the fast variant).
 *   QWE_*   quadrature with extrapolation - always adaptive.
 *   DE_*    double-exponential.  This family contains BOTH a fixed-node method
 *           (Ogata) and an adaptive one (Ooura), and that distinction governs
 *           how the strategy can be used.
 * ------------------------------------------------------------------------- */

/* The names accepted by hankel_transform() live in libhankel.h as
 * LIBHANKEL_ALL_STRATEGIES, so the error message below and the Python one in
 * py_interface.c are built from a single literal and cannot drift apart. */

/* Strategy names for fixed nodes strategies. */
static const char *const FIXED_STRATEGIES =
    "'DHT_Guptasarma', 'DHT_Guptasarma_Fast', 'DHT_Key_51', 'DHT_Key_101', "
    "'DHT_Key_201', 'DHT_Anderson_801', 'Fixed_DE_Ogata'";

/* Map a DHT strategy name to the internal filter index used by
 * hankel_transform_DHT().  These indices are the ones SASfit
 * uses.  Returns 0 if the name is not a DHT strategy. */
static int dht_strategy_index(const char *strategy_name) {
    if (strcmp(strategy_name, "DHT_Guptasarma") == 0) {
        return 6;
    }
    if (strcmp(strategy_name, "DHT_Guptasarma_Fast") == 0) {
        return 7;
    }
    if (strcmp(strategy_name, "DHT_Key_51") == 0) {
        return 8;
    }
    if (strcmp(strategy_name, "DHT_Key_101") == 0) {
        return 9;
    }
    if (strcmp(strategy_name, "DHT_Key_201") == 0) {
        return 10;
    }
    if (strcmp(strategy_name, "DHT_Anderson_801") == 0) {
        return 11;
    }
    return 0;
}

int hankel_transform(int nu, form_factor_f f, double *x, size_t len_x, void *f_ctx, double *output,
                     const char *strategy_name, strategy_params strategy_params) {
    int status;

    if (!(nu == 0 || nu == 1)) {
        fprintf(stderr, "nu needs to be 0 or 1 in order to use the selected strategy\n");
        return -1;
    }

    const int dht_index = dht_strategy_index(strategy_name);
    if (dht_index != 0) {
        for (size_t j = 0; j < len_x; j++) {
            status = hankel_transform_DHT(nu, f, x[j], f_ctx, &output[j], dht_index);
            if (status != 0) {
                return status;
            }
        }
        return 0;
    }

    if (strcmp(strategy_name, "Adaptive_DE_Ooura") == 0) {

        status = validate_n_eval(strategy_params);
        if (status != 0) {
            return status;
        }
        status = validate_eps_rel(strategy_params);
        if (status != 0) {
            return status;
        }

        for (size_t j = 0; j < len_x; j++) {
            status = hankel_transform_DE_Ooura(nu, f, x[j], f_ctx, &output[j],
                                               strategy_params.n_eval, strategy_params.eps_rel);

            if (status != 0) {
                return status;
            }
        }

    } else if (strcmp(strategy_name, "Fixed_DE_Ogata") == 0) {

        status = validate_n_eval(strategy_params);
        if (status != 0) {
            return status;
        }
        status = validate_f_max(strategy_params);
        if (status != 0) {
            return status;
        }

        for (size_t j = 0; j < len_x; j++) {
            status = hankel_transform_DE_Ogata(nu, f, x[j], f_ctx, &output[j],
                                               strategy_params.n_eval, strategy_params.f_max);

            if (status != 0) {
                return status;
            }
        }

    } else if (strcmp(strategy_name, "QWE_Chave") == 0) {

        status = validate_n_eval(strategy_params);
        if (status != 0) {
            return status;
        }
        status = validate_eps_rel(strategy_params);
        if (status != 0) {
            return status;
        }

        for (size_t j = 0; j < len_x; j++) {
            status = hankel_transform_QWE_Chave(nu, f, x[j], f_ctx, &output[j],
                                                strategy_params.n_eval, strategy_params.eps_rel);

            if (status != 0) {
                return status;
            }
        }

    } else if (strcmp(strategy_name, "QWE_Key") == 0) {

        status = validate_n_eval(strategy_params);
        if (status != 0) {
            return status;
        }
        status = validate_eps_rel(strategy_params);
        if (status != 0) {
            return status;
        }

        for (size_t j = 0; j < len_x; j++) {
            status = hankel_transform_QWE_Key(nu, f, x[j], f_ctx, &output[j],
                                              strategy_params.n_eval, strategy_params.eps_rel);

            if (status != 0) {
                return status;
            }
        }

    } else {
        fprintf(stderr, "Invalid strategy name '%s', must be one of : %s.\n", strategy_name,
                LIBHANKEL_ALL_STRATEGIES);
        return -11;
    }

    return 0;
}
