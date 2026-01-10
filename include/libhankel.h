
#ifndef HANKELCLIB_H
#define HANKELCLIB_H

/**
 * Prototypes
 */
// Declare sasfit_params struct
typedef struct
{
	// double             p[50];         
	void *             fparams; 
	double             other_inputs[50];    
	double (*function)(double,  double (*)[50]); //(double, void *);
} hankel_inputs;


/**
 * Function declarations
 */
// // Function pointer type: takes (double, double (*)[50]) -> double
// typedef double (*fn_t)(double, double (*)[50]);
// // Pointer-to-array type
// typedef double (*arr50_t)[50];


// corresponding to strategies 6-11 (NO PARAMS)
double hankel_transform_no_params(
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
	double N_ogata, 
	double h_ogata
);

// Needed for sasfit strategy 0 (DE QUADRATURE)
void sasfit_intdeini(int lenaw, double tiny, double eps, double *aw);
void sasfit_intde(double (*f)(double, void *), double a, double b, double *aw, 
    double *i, double *err, void *fparams);
void sasfit_intdeoini(int lenaw, double tiny, double eps, double *aw);
void sasfit_intdeo(double (*f)(double, void *), double a, double omega, double *aw, 
    double *i, double *err, void *fparams);

double hankel_transform_DE_Quadrature(
	int nu, 
	double (*f)(double, double (*)[50]), 
	double x, 
	double (*fparams)[50], 
	double N_ogata, 
	double eps_nriq
);

// Needed for sasfit strategy 1 (DE QUADRATURE)
double hankel_transform_DE_Ogata(
	int nu, 
	double (*f)(double, double (*)[50]), 
	double x, 
	double (*fparams)[50], 
	double N_ogata, 
	double h_ogata
);


double hankel_transform_QWE_Key(
	int nu, 
	double (*f)(double, double (*)[50]),
	double x, 
	double (*fparams)[50], 
	double N_ogata, 
	double eps_nriq);

double hankel_transform_QWE_Chave(
	int nu, 
	double (*f)(double, double (*)[50]), 
	double x, 
	double (*fparams)[50], 
	double N_ogata, 
	double eps_nriq);


double sasfit_integrate_ctm(
    double int_start,
    double int_end,
    double  (intKern_fct) (double, hankel_inputs *),
    hankel_inputs *param,
    int limit,
    double epsabs,
    double epsrel);

double sasfit_qwe(
	double nu, 
	double (*f)(double, double (*)[50]), 
	double x, 
	void *fparams, 
	int nIntervalsMax, 
	double rtol, 
	double atol
); 
	
double sasfit_HankelChave(
	double order, 
	double (*f)(double, double (*)[50]), 
	double r, 
	void *fparams, 
	int nIntervalsMax, 
	double rerr, 
	double aerr
);

double TanhSinhQuad(
	double (*f)(double, void *), 
	void *par, 
	double a, 
	double b, 
	int n, 
	double eps, 
	double *err
);

double FrJnu(double r, hankel_inputs * inputs);

// function declaration for form factor functions
double form_factor_g_dab(double q, double (*params)[50]); 
double form_factor_sphere(double q, double (*params)[50]);
#endif // HANKELCIB_H

