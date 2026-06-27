
#ifndef LIBHANKEL_H
#define LIBHANKEL_H
#include <stddef.h>

typedef struct {
    int n_eval;      // integer indicating number of function evaluations (N_ogata in SASfit)
    double eps_rel;     // relative error allowed e.g. 1e-9 (eps_nriq in SASfit)
    double f_max;       // float indicating starting guess for max in form factor (h_ogata in SASfit)
} strategy_params;

typedef double (*form_factor_f)(double x, void *ctx);

typedef struct {
    double *params;
} form_factor_ctx;

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
