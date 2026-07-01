
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


int hankel_transform(
    int nu, 
    form_factor_f f, 
    double *x,
    size_t len_x,
    void *f_ctx,
    double *output,
    const char *strategy_name,
    strategy_params strategy_params
); 


// strategies 6-11 (DHT)
double hankel_transform_DHT(
    int nu, 
    form_factor_f f, 
    const double x,
    void *f_ctx,
    double * output,
    int n_strategy
);



// strategy 0 (DE QUADRATURE)
double hankel_transform_DE_Ooura(
    int nu, 
    form_factor_f f, 
    const double x,
    void *f_ctx,
    double * output,
    int n_eval, 
    double eps_rel
);

// strategy 1 (DE QUADRATURE)
double hankel_transform_DE_Ogata(
    int nu, 
    form_factor_f f, 
    const double x,
    void *f_ctx,
    double * output,
    int n_eval, 
    double f_max
);

// strategy 12
double hankel_transform_QWE_Key(
    int nu, 
    form_factor_f f, 
    const double x,
    void *f_ctx,
    double * output,
    int n_eval, 
    double eps_rel
);

// strategy 13
double hankel_transform_QWE_Chave(
    int nu, 
    form_factor_f f, 
    const double x,
    void *f_ctx,
    double * output,
    int n_eval, 
    double eps_rel
);

#endif // LIBHANKEL_H
