#include "libhankel.h"
#include <string.h>



/** 
 * @brief Computes Hankel transform using the method specified by the user.
 * @note The parameter strategy_name can be any of "QWE_Chave", "QWE_Key", 
 *       "DHT_6", "DHT_7", "DHT_8", "DHT_9", "DHT_10", "DHT_11".
 * 
 * @param nu               order of bessel function - must be 0 or 1
 * @param f                pointer to form factor function
 * @param x                value at which to compute the transform
 * @param f_params          params for form factor
 * @param output           pointer to var containing output from transform 
 * @param strategy_name    str corresponding to the strategy name
 * @param strategy_params  struct containing params for the specific strategy
 *                         (please use the docs to find params required by each strategy)
 */ 
double hankel_transform(
    int nu, 
    double (*f)(double, double (*)[50]), 
    double x, 
    double (*f_params)[50], 
    double *output,
    const char *strategy_name,
    strategy_params strategy_params) 
{
    int status;

    if (strcmp(strategy_name, "QWE_Chave") == 0) {
        status = hankel_transform_QWE_Chave(
            nu, 
            f, 
            x, 
            f_params,
            output, 
            strategy_params.n_eval, 
            strategy_params.eps_rel
        );

    } else if (strcmp(strategy_name, "QWE_Key") == 0) {
        status = hankel_transform_QWE_Key(
            nu, 
            f, 
            x, 
            f_params, 
            output,
            strategy_params.n_eval, 
            strategy_params.eps_rel
        );

    } else if (strcmp(strategy_name, "DHT_6") == 0) {
        status = hankel_transform_DHT(
            nu, 
            f, 
            x, 
            f_params, 
            output,
            6
        );

    } else if (strcmp(strategy_name, "DHT_7") == 0) {
        status = hankel_transform_DHT(
            nu, 
            f, 
            x, 
            f_params, 
            output,
            7
        );

    } else if (strcmp(strategy_name, "DHT_8") == 0) {
        status = hankel_transform_DHT(
            nu, 
            f, 
            x, 
            f_params, 
            output,
            8
        );

    } else if (strcmp(strategy_name, "DHT_9") == 0) {
        status = hankel_transform_DHT(
            nu, 
            f, 
            x, 
            f_params, 
            output,
            9
        );

    } else if (strcmp(strategy_name, "DHT_10") == 0) {
        status = hankel_transform_DHT(
            nu, 
            f, 
            x, 
            f_params, 
            output,
            10
        );

    } else if (strcmp(strategy_name, "DHT_11") == 0) {
        status = hankel_transform_DHT(
            nu, 
            f, 
            x, 
            f_params, 
            output,
            11
        );
    }

    return status;
    
}