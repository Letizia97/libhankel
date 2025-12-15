
#include "libhankel.h"

#include "strateg6_const.h"
#include "strateg7_const.h"
#include "strateg8_const.h"
#include "strateg9_const.h"
#include "strateg10_const.h"
#include "strateg11_const.h"

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


double hankel_transform_FBT(double nu, double (*f)(double, double (*)[50]), double x, double (*fparams)[50], double n_ogata, double h_ogata, int n_strategy){
    /* 
    Computes Hankel transform, using FBT.
    Requires additional params - say which - and performance is highly dependent on params.
    Requires knowledge of the form factor function to properly set params. 
    Receives:
        nu         order of bessel function
        *f         pointer to form factor function
        x          value at which to compute the transform
        *fparams   params for form factor
        n_ogata    ... ?
        h_ogata    ... ?
        n_strategy FBT strategy to use --- 1 to 3?
    */

};



double hankel_transform_no_params(double nu, double (*f)(double, double (*)[50]), double x, double (*fparams)[50], int n_strategy){
    /*
    Computes Hankel transform, and does not require additional params.
    Tend to perform well for simple form factors, but it does struggle with oscillatory ones.

    Receives:
        nu         order of bessel function
        *f         pointer to form factor function
        x          value at which to compute the transform
        *fparams   params for form factor
        n_strategy an integer number between 6 and 11
    */
    double res;
    double lambda;
    unsigned int i;
    unsigned int ind;

    ind = nu + 1;

    if (!(nu==0 || nu==1)) {
        fprintf(stderr, "nu needs to be 0 or 1 in order to use the selected strategy\n");
        return 1;
    }

    switch(n_strategy){

        case 6:{
            // HANKEL_GUPTASARMA_97
            res = 0;
            if (nu==0) {
                for (i=0; i<120; i++) {
                    lambda = pow(10.0E0, (aJ0 + i*sJ0)) / x;
                    res = res + (*f)(lambda, fparams) * lambda * WJ0[i] / x;
                }
            } else {
                for (i=0; i<140; i++) {
                    lambda = pow(10.0E0, (aJ1 + i*sJ1)) / x;
                    res = res + (*f)(lambda, fparams) * lambda * WJ1[i] / x;
                }
            }
            break;
        }
        case 7:{
            // HANKEL_GUPTASARMA_97_FAST
            if (nu==0) {
                for (i=0; i<61; i++) {
                    lambda = pow(10.0E0, (aJ0Fast + i*sJ0Fast))/x;
                    res = res+(*f)(lambda, fparams) * lambda * WJ0Fast[i] / x;
                }
            } else {
                for (i=0; i<47; i++) {
                    lambda = pow(10.0E0, (aJ1Fast + i*sJ1Fast)) / x;
                    res = res + (*f)(lambda, fparams) * lambda * WJ0Fast[i] / x;
                }
            }
            break;
        }
        case 8:{
            // HANKEL_KEY_51
            res = 0;
            for (i=0; i<51; i++) {
                lambda = KK51Hankel[i][0] / x;
                res = res + (*f)(lambda, fparams) * lambda * KK51Hankel[i][ind] / x;
            }
            break;
        }
        case 9:{
            // HANKEL_KEY_101
            res = 0;
            for (i=0; i<101; i++) {
                lambda = KK101Hankel[i][0]/x;
                res = res + (*f)(lambda, fparams) * lambda * KK101Hankel[i][ind] / x;
            }
            break;
        }
        case 10:{
            // HANKEL_KEY_201
            res = 0;
            for (i=0; i<201; i++) {
                lambda = KK201Hankel[i][0] / x;
                res = res + (*f)(lambda, fparams) * lambda * KK201Hankel[i][ind] / x;
            }
            break;
        }
        case 11:{
            // HANKEL_ANDERSON_801
            res = 0;
            for (i=0; i<801; i++) {
                lambda = WA801Hankel[i][0] / x;
                res = res + (*f)(lambda, fparams) * lambda * WA801Hankel[i][ind] / x;
            } 
            break;

        }
        default:{
            printf("Strategy number must be integer between 6 and 11\n");
            break;
        }

    }
    return res;
};





// double hankel_strategy_6(double nu, double (*f)(double, double (*)[50]), double x, double (*fparams)[50]){
//     /*
//     Computes Hankel transform, strategy 7 (HANKEL_GUPTASARMA_97 in SASfit)
//     */
//     double res, lambda;
//     unsigned int i;

//     if (!(nu==0 || nu==1)) {
//         fprintf(stderr, "nu needs to be 0 or 1 in order to use hankel_strategy_6\n");
//         return 1;
//     }

//     if (nu==0) {
//         res = 0;
//         for (i=0; i<120; i++) {
//             lambda = pow(10.0E0,(aJ0 + i*sJ0)) / x;
//             res = res + (*f)(lambda,fparams) * lambda * WJ0[i] / x;
//         }
//     } else {
//         res = 0;
//         for (i=0; i<140; i++) {
//             lambda = pow(10.0E0,(aJ1 + i*sJ1)) / x;
//             res = res + (*f)(lambda,fparams) * lambda * WJ1[i] / x;
//         }
//     }
//     return res;
// };


// double hankel_strategy_8(double nu, double (*f)(double, double (*)[50]), double x, double (*fparams)[50]){
//     /*
//     Computes Hankel transform, strategy 8 (HANKEL_KEY_51 in SASfit)
//     */
//     double res, lambda;
//     unsigned int i;

//     if (!(nu==0 || nu==1)) {
//         fprintf(stderr, "nu needs to be 0 or 1 in order to use hankel_strategy_8\n");
//         return 1;
//     }
//     res = 0;
//     lambda = KK51Hankel[i][0] / x;
//     if (nu==0) {
        
//         for (i=0; i<51; i++) {
//             res = res + (*f)(lambda,fparams) * lambda * KK51Hankel[i][1] / x;
//         }
//     } else {
//         for (i=0; i<51; i++) {
//             res = res + (*f)(lambda,fparams) * lambda * KK51Hankel[i][2] / x;
//         }
//     }
//     return res;
  
// };



// double hankel_strategy_9(double nu, double (*f)(double, double (*)[50]), double x, double (*fparams)[50]){
//     /*
//     Computes Hankel transform, strategy 9 (HANKEL_KEY_101 in SASfit)
//     */
//     double res, lambda;
//     unsigned int i;

//     if (!(nu==0 || nu==1)) {
//         fprintf(stderr, "nu needs to be 0 or 1 in order to use hankel_strategy_9\n");
//         return 1;
//     }
//     if (nu==0) {
//         res = 0;
//         for (i=0; i<101; i++) {
//             lambda = KK101Hankel[i][0]/x;
//             res = res + (*f)(lambda,fparams) * lambda * KK101Hankel[i][1] / x;
//         }
//     } else {
//         res = 0;
//         for (i=0; i<101; i++) {
//             lambda = KK101Hankel[i][0]/x;
//             res = res + (*f)(lambda,fparams) * lambda * KK101Hankel[i][2] / x;
//         }
//     }

//     return res;
  
// };


// double hankel_strategy_10(double nu, double (*f)(double, double (*)[50]), double x, double (*fparams)[50]){
//     /*
//     Computes Hankel transform, strategy 10 (HANKEL_KEY_201 in SASfit)

//     Receives:
//         nu       double
//         f        function with signature: double (*f)(double, void *)
//         x        double, value where to compute the transform
//         fparams  void pointer to parameters to be passed to f
    
//     Returns:
//         res      
    
//     */
//     double res, lambda;
//     unsigned int i;

//     if (!(nu==0 || nu==1)) {
//         fprintf(stderr, "nu needs to be 0 or 1 in order to use hankel_strategy_10\n");
//         return 1;
//     }
//     if (nu==0) {
//         res = 0;
//             for (i=0; i<201; i++) {
//                 lambda = KK201Hankel[i][0] / x;
//                 res = res + (*f)(lambda,fparams) * lambda * KK201Hankel[i][1] / x;
//             }
//     } else {
//         res = 0;
//             for (i=0; i<201; i++) {
//                 lambda = KK201Hankel[i][0] / x;
//                 res = res + (*f)(lambda,fparams) * lambda * KK201Hankel[i][2] / x;
//             }
//     }

//     return res;
// };


