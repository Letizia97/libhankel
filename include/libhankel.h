
#ifndef HANKELCLIB_H
#define HANKELCLIB_H

/**
 * Prototypes
 */
// Declare sasfit_params struct
// typedef struct
// {
// 	double             p[MAXPAR];         //!< Parameter of a function.
// 	double             *xarr;             //!< hack for OZ solver
// 	double             *yarr;             //!< hack for OZ solver
// 	void *             moreparam; 
// 	double             more_p[MAXPAR];    //!< more Parameter of a function.
// 	double (*function)(double, void *);
// } sasfit_param;


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

// function declaration for form factor functions
double form_factor_g_dab(double q, double (*params)[50]); 
double form_factor_sphere(double q, double (*params)[50]);
#endif // HANKELCIB_H

