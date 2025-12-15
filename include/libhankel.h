
#ifndef HANKELCLIB_H
#define HANKELCLIB_H


/**
 * Constants
 */
// #define MAXPAR 50

// Declarations: visible to other files
// extern const double aJ0, sJ0, aJ1, sJ1;
// extern const double aJ1Fast, sJ1Fast, aJ0Fast, sJ0Fast;
// extern const double WJ0[120];
// extern const double WJ1[140];
// extern const double WJ0Fast[61];
// extern const double WJ1Fast[47];
// extern const double KK51Hankel[51][3];
// extern const double KK101Hankel[101][3];
// extern const double KK201Hankel[201][3];

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

// function declaration for form factor functions
double sasfit_ff_g_dab(double q, double (*params)[50]);
double sasfit_ff_sphere(double q, double (*params)[50]);

// // function declaration for hankel tr, strategy 9
// double hankel_strategy_7(
// 	double nu, 
// 	double (*f)(double, double (*)[50]), 
// 	double x, 
// 	double (*fparams)[50]
// );

// // function declaration for hankel tr, strategy 9
// double hankel_strategy_8(
// 	double nu, 
// 	double (*f)(double, double (*)[50]), 
// 	double x, 
// 	double (*fparams)[50]
// );

// // function declaration for hankel tr, strategy 9
// double hankel_strategy_9(
// 	double nu, 
// 	double (*f)(double, double (*)[50]), 
// 	double x, 
// 	double (*fparams)[50]
// );

// // function declaration for hankel tr, strategy 10
// double hankel_strategy_10(
// 	double nu, 
// 	double (*f)(double, double (*)[50]), 
// 	double x, 
// 	double (*fparams)[50]
// );


double hankel_transform_no_params(
	double nu, 
	double (*f)(double, double (*)[50]), 
	double x, 
	double (*fparams)[50], 
	int n_strategy
);

#endif // HANKELCIB_H

