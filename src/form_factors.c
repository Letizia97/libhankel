#include "include/libhankel.h"
#include <math.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_bessel.h>

double form_factor_g_dab(double q, double (*params)[50]) {
    /*
    Compute the g_dab form factor. 

    Receives: 
        - double q          value at which to compute the Hankel t
        - double (*params)  pointer to an array of params for the function

    */
    double XI = (*params)[0];
    double H = (*params)[1];
    double ETA = (*params)[2];
    double numer, denom;

    numer = gsl_pow_2(gsl_pow_3(2*XI)*gsl_sf_poch(H,1.5)*ETA)*gsl_pow_3(M_PI);
    denom = pow(1+gsl_pow_2(q*XI),1.5+H);
    return numer / denom;
};


double form_factor_sphere(double q, double (*params)[50]) {
    /*
    Compute the sphere form factor. 

    Receives: 
        - double q          value at which to compute the Hankel t
        - double (*params)  pointer to an array of params for the function

    */
    double R = (*params)[0];
    double ETA = (*params)[1];
    double interm;

	if (q * R < 1e-4) {
		interm = ETA*4.0/3.0*M_PI*R*R*R*(1 - gsl_pow_2(q*R)/10. + gsl_pow_4(q*R)/280. - gsl_pow_6(q*R)/15120.);
	} else {
		interm = ETA*4.0*M_PI*(sin(q*R) - q*R*cos(q*R))/gsl_pow_3(q);
	}
    return gsl_pow_2(interm);

};