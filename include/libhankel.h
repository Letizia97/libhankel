
#ifndef LIBHANKEL_H
#define LIBHANKEL_H


typedef struct {
    double n_eval;      // integer indicating number of function evaluations (N_ogata in SASfit)
    double eps_rel;     // relative error allowed e.g. 1e-9 (eps_nriq in SASfit)
    double f_max;       // float indicating starting guess for max in form factor (h_ogata in SASfit)
} strategy_params;


double hankel_transform(
    int nu, 
    double (*f)(double, double (*)[50]), 
    double x, 
    double (*f_params)[50], 
    double *output,
    const char *strategy_name,
    strategy_params strategy_params
); 


// strategies 6-11 (DHT)
double hankel_transform_DHT(
	int nu, 
	double (*f)(double, double (*)[50]), 
	double x, 
	double (*f_params)[50], 
	double (*output),
	int n_strategy
);

// corresponding to strategies 2-4 (FBT)
double compute_hankel_FBT(
	int nu, 
	double (*f)(double, double (*)[50]), 
	double x, 
	double (*f_params)[50], 
	int n_method,
	int N, 
	double h
);
double hankel_transform_FBT(
	int nu, 
	double (*f)(double, double (*)[50]), 
	double x, 
	double (*f_params)[50], 
	int n_method, 
	double n_eval, 
	double f_max
);

// strategy 0 (DE QUADRATURE)
double hankel_transform_DE_Quadrature(
    int nu, 
    double (*f)(double, double (*)[50]), 
    double x, 
    double (*f_params)[50], 
    double *output,
    double n_eval, 
    double eps_rel
);

// strategy 1 (DE QUADRATURE)
double hankel_transform_DE_Ogata(
    int nu, 
    double (*f)(double, double (*)[50]), 
    double x, 
    double (*f_params)[50], 
    double *output,
    double n_eval, 
    double f_max
);

// strategy 12
double hankel_transform_QWE_Key(
	int nu, 
	double (*f)(double, double (*)[50]),
	double x, 
	double (*f_params)[50], 
	double (*output),
	double n_eval, 
	double eps_rel
);

// strategy 13
double hankel_transform_QWE_Chave(
	int nu, 
	double (*f)(double, double (*)[50]), 
	double x, 
	double (*f_params)[50], 
	double (*output),
	double n_eval, 
	double eps_rel
);

#endif // LIBHANKEL_H
