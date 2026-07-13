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

int hankel_transform(int nu, form_factor_f f, double *x, size_t len_x, void *f_ctx, double *output,
                     const char *strategy_name, strategy_params strategy_params) {
    int status;

    if (!(nu == 0 || nu == 1)) {
        fprintf(stderr, "nu needs to be 0 or 1 in order to use the selected strategy\n");
        return -1;
    }

    if (strcmp(strategy_name, "DE_Ooura") == 0) {

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

    } else if (strcmp(strategy_name, "DE_Ogata") == 0) {

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

    } else if (strcmp(strategy_name, "DHT_6") == 0) {

        for (size_t j = 0; j < len_x; j++) {
            status = hankel_transform_DHT(nu, f, x[j], f_ctx, &output[j], 6);

            if (status != 0) {
                return status;
            }
        }

    } else if (strcmp(strategy_name, "DHT_7") == 0) {
        for (size_t j = 0; j < len_x; j++) {
            status = hankel_transform_DHT(nu, f, x[j], f_ctx, &output[j], 7);

            if (status != 0) {
                return status;
            }
        }

    } else if (strcmp(strategy_name, "DHT_8") == 0) {
        for (size_t j = 0; j < len_x; j++) {
            status = hankel_transform_DHT(nu, f, x[j], f_ctx, &output[j], 8);
            if (status != 0) {
                return status;
            }
        }

    } else if (strcmp(strategy_name, "DHT_9") == 0) {
        for (size_t j = 0; j < len_x; j++) {
            status = hankel_transform_DHT(nu, f, x[j], f_ctx, &output[j], 9);
            if (status != 0) {
                return status;
            }
        }

    } else if (strcmp(strategy_name, "DHT_10") == 0) {
        for (size_t j = 0; j < len_x; j++) {
            status = hankel_transform_DHT(nu, f, x[j], f_ctx, &output[j], 10);
            if (status != 0) {
                return status;
            }
        }

    } else if (strcmp(strategy_name, "DHT_11") == 0) {
        for (size_t j = 0; j < len_x; j++) {
            status = hankel_transform_DHT(nu, f, x[j], f_ctx, &output[j], 11);
            if (status != 0) {
                return status;
            }
        }
    } else {
        fprintf(stderr, "strategy_name must be one of the strings specified in the docs\n");
        return -11;
    }
}