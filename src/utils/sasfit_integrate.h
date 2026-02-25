
#ifndef UTILS_H
#define UTILS_H


typedef struct
{
	void *             f_params; 
	double             other_inputs[50];    
	double (*function)(double,  double (*)[50]); //(double, void *);
} hankel_inputs;

// Needed for sasfit strategy 12 and 13
double sasfit_integrate_ctm(
    double int_start,
    double int_end,
    double  (intKern_fct) (double, hankel_inputs *),
    hankel_inputs *param,
    int limit,
    double epsabs,
    double epsrel);


// input to sasfit_integrate_ctm
double FrJnu(double r, hankel_inputs * inputs);


#endif // UTILS_H

