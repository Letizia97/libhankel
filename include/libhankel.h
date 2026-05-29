
#ifndef LIBHANKEL_H
#define LIBHANKEL_H
#include <stddef.h>

typedef struct {
    int n_eval;      // integer indicating number of function evaluations (N_ogata in SASfit)
    double eps_rel;     // relative error allowed e.g. 1e-9 (eps_nriq in SASfit)
    double f_max;       // float indicating starting guess for max in form factor (h_ogata in SASfit)
} strategy_params;

typedef double (*intern_fct_type)(double x, void *ctx);

typedef double (*form_factor_f)(double x, double *params, size_t n_params, void *user_data);

int hankel_transform(
    int nu, 
    form_factor_f f, 
    double *x,
    size_t len_x,
    double *f_params,
    size_t n_params,
    double * output,
    const char *strategy_name,
    strategy_params strategy_params,
    void *user_data
); 


// strategies 6-11 (DHT)
double hankel_transform_DHT(
    int nu, 
    form_factor_f f, 
    const double x,
    double *f_params,
    size_t n_params,
    double * output,
    int n_strategy,
    void *user_data
);

// corresponding to strategies 2-4 (FBT)
double compute_hankel_FBT(
    int nu, 
    intern_fct_type intern_fct,
    void *ctx,
    double x, 
    int n_method, 
    int n_eval, 
    double f_max,
    void *user_data
);
double hankel_transform_FBT(
    int nu, 
    double (*f)(double, double (*)[50]), 
    double x, 
    double (*f_params)[50], 
    double *output,
    int n_method, 
    int n_eval, 
    double f_max,
    void *user_data
);

// strategy 0 (DE QUADRATURE)
double hankel_transform_DE_Ooura(
    int nu, 
    form_factor_f f, 
    const double x,
    double *f_params,
    size_t n_params,
    double * output,
    int n_eval, 
    double eps_rel,
    void *user_data
);

// strategy 1 (DE QUADRATURE)
double hankel_transform_DE_Ogata(
    int nu, 
    form_factor_f f, 
    const double x,
    double *f_params,
    size_t n_params,
    double * output,
    int n_eval, 
    double f_max,
    void *user_data
);

// strategy 12
double hankel_transform_QWE_Key(
    int nu, 
    form_factor_f f, 
    const double x,
    double *f_params,
    size_t n_params,
    double * output,
    int n_eval, 
    double eps_rel,
    void *user_data
);

// strategy 13
double hankel_transform_QWE_Chave(
    int nu, 
    form_factor_f f, 
    const double x,
    double *f_params,
    size_t n_params,
    double * output,
    int n_eval, 
    double eps_rel,
    void *user_data
);

#endif // LIBHANKEL_H
