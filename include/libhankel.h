
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
     *        Corresponds to ``N_ogata`` in SASfit.
     */
    int n_eval;

    /**
     * @brief Tolerance parameter or relative error allowed.
     *
     * Determines the numerical precision or convergence threshold
     * of the computation. Smaller values lead to higher accuracy
     * at the cost of increased computation time.
     * Corresponds to ``eps_nriq`` in SASfit.
     */
    double eps_rel;

    /**
     * @brief Starting guess for the maximum in the function \f$ x \cdot formFactor \f$.
     *        Corresponds to ``h_ogata`` in SASfit.
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
 * @def LIBHANKEL_ALL_STRATEGIES
 * @brief Every strategy name accepted by @ref hankel_transform, as a quoted,
 *        comma-separated list.
 *
 * Defined here so the C and the Python error messages are built from the same
 * literal and cannot drift apart when a strategy is added or renamed.
 */
#define LIBHANKEL_ALL_STRATEGIES                                                                   \
    "'DHT_Guptasarma', 'DHT_Guptasarma_Fast', 'DHT_Key_51', 'DHT_Key_101', "                       \
    "'DHT_Key_201', 'DHT_Anderson_801', 'Fixed_DE_Ogata', 'Adaptive_DE_Ooura', "                   \
    "'QWE_Chave', 'QWE_Key'"

/**
 * @brief Computes Hankel transform using the method specified by the user.
 *
 * @note strategy_name can be any of
 *       "DHT_Guptasarma", "DHT_Guptasarma_Fast", "DHT_Key_51", "DHT_Key_101",
 *       "DHT_Key_201", "DHT_Anderson_801", "Fixed_DE_Ogata", "Adaptive_DE_Ooura",
 *       "QWE_Chave", "QWE_Key" (see @ref LIBHANKEL_ALL_STRATEGIES).
 *
 * @note The digital filters are named after the filter they implement -
         the three authors are independent designs, and accuracy is not
 *       monotonic in the tap count.
 *
 * @param nu               order of bessel function - must be 0 or 1
 * @param f                pointer to function to transform (see @ref form_factor_f).
 * @param x                pointer to array of x at which to compute the transform
 * @param f_ctx            pointer to struct containing inputs for f
 * @param output           pointer to array containing output from transform
 * @param strategy_name    str corresponding to the strategy name
 * @param strategy_params  struct containing params for the specific strategy.
 *                         See @ref strategy_params for all elements this struct
 *                         can contain and the page
 *                         <a href="../usage/strategy_params.html">Strategy Parameters</a>
 *                         to check the params required by each strategy.
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
 *
 * @return 0 on success, or a negative status code (see the
 *         <a href="../usage/status_codes.html">Status Codes</a> page). @p x
 *         must be finite and greater than zero.
 */
int hankel_transform_DHT(int nu, form_factor_f f, double x, void *f_ctx, double *output,
                         int n_strategy);

/**
 * @brief Computes Hankel transform, using de-quadrature.
 * @note Corresponds to strategy 0 in SASfit, or HANKEL_OOURA_DEO.
 *
 * Takes the whole array at once: the node/weight tables built by Ooura's
 * initialisers dominate the cost and depend only on @p eps_rel, so one call
 * builds them once instead of once per point.
 *
 * @param nu         order of bessel function - must be 0 or 1
 * @param f          pointer to function to transform
 * @param x          pointer to array of values at which to compute the transform
 * @param len_x      number of entries in @p x and in @p output
 * @param f_ctx      pointer to struct containing inputs for f
 * @param output     pointer to array containing output from transform
 * @param n_eval     integer indicating number of function evaluations (``N_ogata`` in SASfit)
 * @param eps_rel    relative error allowed e.g. 1e-9 (``eps_nriq`` in SASfit)
 *
 * @return 0 on success, or a negative status code (see the
 *         <a href="../usage/status_codes.html">Status Codes</a> page). Every
 *         entry of @p x must be finite and greater than zero; the whole array
 *         is checked before any work is done.
 */
int hankel_transform_DE_Ooura(int nu, form_factor_f f, const double *x, size_t len_x, void *f_ctx,
                              double *output, int n_eval, double eps_rel);

/**
 * @brief Computes Hankel transform, using de-quadrature.
 * @note Corresponds to strategy 1 in SASfit or HANKEL_OGATA_2005.
 *
 * Takes the whole array at once: everything but the form factor evaluation
 * depends only on @p nu, the node index and @p f_max - including the Bessel
 * zeros, whose root search dominates the runtime - so one call builds the node
 * table once instead of once per point.
 *
 * @param nu         order of bessel function - must be 0 or 1
 * @param f          pointer to function to transform
 * @param x          pointer to array of values at which to compute the transform
 * @param len_x      number of entries in @p x and in @p output
 * @param f_ctx      pointer to struct containing inputs for f
 * @param output     pointer to array containing output from transform
 * @param n_eval     integer indicating number of function evaluations (``N_ogata`` in SASfit)
 * @param f_max      float indicating starting guess for max in form factor (``h_ogata`` in SASfit)
 *
 * @return 0 on success, or a negative status code (see the
 *         <a href="../usage/status_codes.html">Status Codes</a> page). Every
 *         entry of @p x must be finite and greater than zero; the whole array
 *         is checked before any work is done.
 */
int hankel_transform_DE_Ogata(int nu, form_factor_f f, const double *x, size_t len_x, void *f_ctx,
                              double *output, int n_eval, double f_max);

/**
 * @brief Computes Hankel transform using the Quadrature With Extrapolation method by Key.
 * @note Corresponds to strategy 12 in SASfit.
 *
 * @param nu         order of bessel function - must be 0 or 1
 * @param f          pointer to function to transform
 * @param x          value at which to compute the transform
 * @param f_ctx      pointer to struct containing inputs for f
 * @param output     pointer to var containing output from transform
 * @param n_eval     integer indicating number of function evaluations (``N_ogata`` in SASfit)
 * @param eps_rel    relative error allowed e.g. 1e-9 (``eps_nriq`` in SASfit)
 *
 * @return 0 on success, or a negative status code (see the
 *         <a href="../usage/status_codes.html">Status Codes</a> page). @p x
 *         must be finite and greater than zero.
 */
int hankel_transform_QWE_Key(int nu, form_factor_f f, double x, void *f_ctx, double *output,
                             int n_eval, double eps_rel);

/**
 * @brief Computes Hankel transform using the Quadrature With Extrapolation method by Chave.
 * @note Corresponds to strategy 13 in SASfit.
 *
 * @param nu         order of bessel function - must be 0 or 1
 * @param f          pointer to function to transform
 * @param x          value at which to compute the transform
 * @param f_ctx      pointer to struct containing inputs for f
 * @param output     pointer to var containing output from transform
 * @param n_eval     integer indicating number of function evaluations (``N_ogata`` in SASfit)
 * @param eps_rel    relative error allowed e.g. 1e-9 (``eps_nriq`` in SASfit)
 *
 * @return 0 on success, or a negative status code (see the
 *         <a href="../usage/status_codes.html">Status Codes</a> page). @p x
 *         must be finite and greater than zero.
 */
int hankel_transform_QWE_Chave(int nu, form_factor_f f, double x, void *f_ctx, double *output,
                               int n_eval, double eps_rel);

#endif // LIBHANKEL_H
