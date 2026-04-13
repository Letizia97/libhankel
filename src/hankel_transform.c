#include "libhankel.h"
#include <string.h>
#include <stddef.h>
#include <stdio.h>


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
        return -8;
    }
    return 0;
}

int validate_f_max(strategy_params strategy_params) {
    if (strategy_params.f_max == 0) {
        fprintf(stderr, "Error: f_max must be provided and cannot be zero\n");
        return -8;
    }
    return 0;
}


/** 
 * @brief Computes Hankel transform using the method specified by the user.
 * @note The parameter strategy_name can be any of "QWE_Chave", "QWE_Key", 
 *       "DHT_6", "DHT_7", "DHT_8", "DHT_9", "DHT_10", "DHT_11".
 * 
 * @param nu               order of bessel function - must be 0 or 1
 * @param f                pointer to form factor function
 * @param x                pointer to array of x values at which to compute the transform
 * @param f_params         params for form factor
 * @param output           pointer to array containing output from transform 
 * @param len_x            length of array x 
 * @param strategy_name    str corresponding to the strategy name
 * @param strategy_params  struct containing params for the specific strategy
 *                         (please use the docs to find params required by each strategy)
 */ 
int hankel_transform(
    int nu, 
    double (*f)(double, double (*)[50]), 
    const double *x,
    double (*f_params)[50], 
    double * output,
    int len_x,
    const char *strategy_name,
    strategy_params strategy_params) 
{
    int status;

    if (!(nu==0 || nu==1)) {
        fprintf(stderr, "nu needs to be 0 or 1 in order to use the selected strategy\n");
        return -1;
    }

    if (strcmp(strategy_name, "DE_Ooura") == 0) {

        status = validate_n_eval(strategy_params);
        if (status!=0) {return status;}
        status = validate_eps_rel(strategy_params);
        if (status!=0) {return status;}

        for (int j = 0; j < len_x; j++) {
            status = hankel_transform_DE_Ooura(
                nu, 
                f, 
                x[j], 
                f_params,
                &output[j], 
                strategy_params.n_eval, 
                strategy_params.eps_rel
            );

            if (status!=0) {return status;}
        }

    } else if (strcmp(strategy_name, "DE_Ogata") == 0) {

        status = validate_n_eval(strategy_params);
        if (status!=0) {return status;}
        status = validate_f_max(strategy_params);
        if (status!=0) {return status;}

        for (int j = 0; j < len_x; j++) {
            status = hankel_transform_DE_Ogata(
                nu, 
                f, 
                x[j], 
                f_params,
                &output[j], 
                strategy_params.n_eval, 
                strategy_params.f_max
            );

            if (status!=0) {return status;}
        }

    } else if (strcmp(strategy_name, "QWE_Chave") == 0) {

        status = validate_n_eval(strategy_params);
        if (status!=0) {return status;}
        status = validate_eps_rel(strategy_params);
        if (status!=0) {return status;}

        for (int j = 0; j < len_x; j++) {
            status = hankel_transform_QWE_Chave(
                nu, 
                f, 
                x[j], 
                f_params,
                &output[j], 
                strategy_params.n_eval, 
                strategy_params.eps_rel
            );

            if (status!=0) {return status;}
        }

    } else if (strcmp(strategy_name, "QWE_Key") == 0) {

        status = validate_n_eval(strategy_params);
        if (status!=0) {return status;}
        status = validate_eps_rel(strategy_params);
        if (status!=0) {return status;}

        for (int j = 0; j < len_x; j++) {
            status = hankel_transform_QWE_Key(
                nu, 
                f, 
                x[j], 
                f_params, 
                &output[j],
                strategy_params.n_eval, 
                strategy_params.eps_rel
            );

            if (status!=0) {return status;}
        }

    } else if (strcmp(strategy_name, "DHT_6") == 0) {

        for (int j = 0; j < len_x; j++) {
            status = hankel_transform_DHT(
                nu, 
                f, 
                x[j], 
                f_params, 
                &output[j],
                6
            );

            if (status!=0) {return status;}
        }

    } else if (strcmp(strategy_name, "DHT_7") == 0) {
        for (int j = 0; j < len_x; j++) {
            status = hankel_transform_DHT(
                nu, 
                f, 
                x[j], 
                f_params, 
                &output[j],
                7
            );

            if (status!=0) {return status;}
        }

    } else if (strcmp(strategy_name, "DHT_8") == 0) {
        for (int j = 0; j < len_x; j++) {
            status = hankel_transform_DHT(
                nu, 
                f, 
                x[j], 
                f_params, 
                &output[j],
                8
            );
            if (status!=0) {return status;}
        }

    } else if (strcmp(strategy_name, "DHT_9") == 0) {
        for (int j = 0; j < len_x; j++) {
            status = hankel_transform_DHT(
                nu, 
                f, 
                x[j], 
                f_params, 
                &output[j],
                9
            );
            if (status!=0) {return status;}
        }

    } else if (strcmp(strategy_name, "DHT_10") == 0) {
        for (int j = 0; j < len_x; j++) {
            status = hankel_transform_DHT(
                nu, 
                f, 
                x[j], 
                f_params, 
                &output[j],
                10
            );
            if (status!=0) {return status;}
        }

    } else if (strcmp(strategy_name, "DHT_11") == 0) {
        for (int j = 0; j < len_x; j++) {
            status = hankel_transform_DHT(
                nu, 
                f, 
                x[j], 
                f_params, 
                &output[j],
                11
            );
            if (status!=0) {return status;}
        }
    }

    return status;
    
}