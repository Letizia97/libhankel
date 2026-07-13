
#ifndef LIBHANKEL_H
#define LIBHANKEL_H
#include <stddef.h>

/**
 * @struct strategy_params
 * @brief Parameters required by the selected Hankel strategy.
 *        This structure groups together parameters that influence
 *        the behavior and accuracy of the algorithm.
 */
typedef struct {

    /**
     * @brief Integer parameter indicating the maximum allowed
     *        number of function evaluations.
     *        Corresponds to N_ogata in SASfit.
     */
    int n_eval;

    /**
     * @brief Tolerance parameter or relative error allowed.
     *
     * Determines the numerical precision or convergence threshold
     * of the computation. Smaller values lead to higher accuracy
     * at the cost of increased computation time.
     * Corresponds to eps_nriq in SASfit.
     */
    double eps_rel;

    /**
     * @brief Starting guess for the maximum in the function x * form_factor.
     *        Corresponds to h_ogata in SASfit.
     */
    double f_max;
} strategy_params;

/**
 * @typedef form_factor_f
 * @brief Function pointer type for a form factor callback.
 *
 * This callback is used to evaluate the form factor at a given input value.
 *
 * @param x   Input variable .
 * @param ctx User-provided context or data pointer (may be NULL).
 *
 * @return The value of the form factor evaluated at @p x.
 */
typedef double (*form_factor_f)(double x, void *ctx);

/**
 * @brief Computes Hankel transform using the method specified by the user.
 * @note The parameter strategy_name can be any of "QWE_Chave", "QWE_Key",
 *       "DHT_6", "DHT_7", "DHT_8", "DHT_9", "DHT_10", "DHT_11".
 *
 * @param nu               order of bessel function - must be 0 or 1
 * @param f                pointer to function to transform (see @ref
 * form_factor_f).
 * @param x                pointer to array of x at which to compute the
 * transform
 * @param f_ctx            pointer to struct containing inputs for f
 * @param output           pointer to array containing output from transform
 * @param strategy_name    str corresponding to the strategy name
 * @param strategy_params  struct containing params for the specific strategy.
 *                         See @ref strategy_params for all elements this struct
 * can contain and the page <a href="../usage/strategy_params.html">Strategy
 * Parameters</a> to check the params required by each strategy.
 */
int hankel_transform(int nu, form_factor_f f, double *x, size_t len_x, void *f_ctx, double *output,
                     const char *strategy_name, strategy_params strategy_params);

/**
 * @brief Computes Hankel transform, using digital filters.
 * @note Corresponds to strategies 6-11 in SASfit.
 * @note Does not allow any error control.
 * @note Perform wells for simple form factors, but it does struggle
 *       with oscillatory ones. For oscillatory form factors, it is
 *       advisable to start with this method for a rough Hankel transform
 *       computation, and then refine it with "hankel_transform_QWE_Chave"
 *       or "hankel_transform_QWE_Key".
 *
 * @param nu         order of bessel function - must be 0 or 1
 * @param f          pointer to function to transform
 * @param x          value at which to compute the transform
 * @param f_ctx      pointer to struct containing inputs for f
 * @param output     pointer to var containing output from transform
 * @param n_strategy an integer number between 6 and 11
 *                   (determines which weightings to use for the transform
 *                   and corresponds to SASfit strategies 6-11)
 */
double hankel_transform_DHT(int nu, form_factor_f f, const double x, void *f_ctx, double *output,
                            int n_strategy);

/**
 * @brief Computes Hankel transform, using de-quadrature.
 * @note Corresponds to strategy 0 in SASfit, or HANKEL_OOURA_DEO.
 *
 * @param nu         order of bessel function - must be 0 or 1
 * @param f          pointer to function to transform
 * @param x          value at which to compute the transform
 * @param f_ctx      pointer to struct containing inputs for f
 * @param output     pointer to var containing output from transform
 * @param n_eval     integer indicating number of function evaluations (N_ogata
 * in SASfit)
 * @param eps_rel    relative error allowed e.g. 1e-9 (eps_nriq in SASfit)
 */
double hankel_transform_DE_Ooura(int nu, form_factor_f f, const double x, void *f_ctx,
                                 double *output, int n_eval, double eps_rel);

/**
 * @brief Computes Hankel transform, using de-quadrature.
 * @note Corresponds to strategy 1 in SASfit or HANKEL_OGATA_2005.
 *
 * @param nu         order of bessel function - must be 0 or 1
 * @param f          pointer to function to transform
 * @param x          value at which to compute the transform
 * @param f_ctx      pointer to struct containing inputs for f
 * @param output     pointer to var containing output from transform
 * @param n_eval     integer indicating number of function evaluations (N_ogata
 * in SASfit)
 * @param f_max      float indicating starting guess for max in form factor
 * (h_ogata in SASfit)
 */
double hankel_transform_DE_Ogata(int nu, form_factor_f f, const double x, void *f_ctx,
                                 double *output, int n_eval, double f_max);

/**
 * @brief Computes Hankel transform using the Quadrature With Extrapolation
 * method by Key.
 * @note Corresponds to strategy 12 in SASfit.
 *
 * @param nu         order of bessel function - must be 0 or 1
 * @param f          pointer to function to transform
 * @param x          value at which to compute the transform
 * @param f_ctx      pointer to struct containing inputs for f
 * @param output     pointer to var containing output from transform
 * @param n_eval     integer indicating number of function evaluations (N_ogata
 * in SASfit)
 * @param eps_rel    relative error allowed e.g. 1e-9 (eps_nriq in SASfit)
 */
double hankel_transform_QWE_Key(int nu, form_factor_f f, const double x, void *f_ctx,
                                double *output, int n_eval, double eps_rel);

/**
 * @brief Computes Hankel transform using the Quadrature With Extrapolation
 * method by Chave.
 * @note Corresponds to strategy 13 in SASfit.
 *
 * @param nu         order of bessel function - must be 0 or 1
 * @param f          pointer to function to transform
 * @param x          value at which to compute the transform
 * @param f_ctx      pointer to struct containing inputs for f
 * @param output     pointer to var containing output from transform
 * @param n_eval     integer indicating number of function evaluations (N_ogata
 * in SASfit)
 * @param eps_rel    relative error allowed e.g. 1e-9 (eps_nriq in SASfit)
 */
double hankel_transform_QWE_Chave(int nu, form_factor_f f, const double x, void *f_ctx,
                                  double *output, int n_eval, double eps_rel);

#endif // LIBHANKEL_H
