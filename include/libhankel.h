
#ifndef LIBHANKEL_H
#define LIBHANKEL_H


// strategies 6-11 (DHT)
double hankel_transform_DHT(
	int nu, 
	double (*f)(double, double (*)[50]), 
	double x, 
	double (*fparams)[50], 
	int n_strategy
);

// corresponding to strategies 2-4 (FBT)
double compute_hankel_FBT(
	int nu, 
	double (*f)(double, double (*)[50]), 
	double x, 
	double (*fparams)[50], 
	int n_method,
	int N, 
	double h
);
double hankel_transform_FBT(
	int nu, 
	double (*f)(double, double (*)[50]), 
	double x, 
	double (*fparams)[50], 
	int n_method, 
	double n_eval, 
	double f_max
);

// strategy 0 (DE QUADRATURE)
double hankel_transform_DE_Quadrature(
	int nu, 
	double (*f)(double, double (*)[50]), 
	double x, 
	double (*fparams)[50], 
	double n_eval, 
	double eps_rel
);

// strategy 1 (DE QUADRATURE)
double hankel_transform_DE_Ogata(
	int nu, 
	double (*f)(double, double (*)[50]), 
	double x, 
	double (*fparams)[50], 
	double n_eval, 
	double f_max
);

// strategy 12
double hankel_transform_QWE_Key(
	int nu, 
	double (*f)(double, double (*)[50]),
	double x, 
	double (*fparams)[50], 
	double n_eval, 
	double eps_rel
);

// strategy 13
double hankel_transform_QWE_Chave(
	int nu, 
	double (*f)(double, double (*)[50]), 
	double x, 
	double (*fparams)[50], 
	double n_eval, 
	double eps_rel
);

// function declaration for form factor functions
double form_factor_g_dab(double q, double (*params)[50]); 
double form_factor_sphere(double q, double (*params)[50]);
double form_factor_broad_peak(double q, double (*params)[50]);

#endif // LIBHANKEL_H

