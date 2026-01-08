
#include "include/libhankel.h"
#include "include/strateg6_const.h"
#include "include/strateg7_const.h"
#include "include/strateg8_const.h"
#include "include/strateg9_const.h"
#include "include/strateg10_const.h"
#include "include/strateg11_const.h"

#include <math.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_bessel.h>


double hankel_transform_no_params(int nu, double (*f)(double, double (*)[50]), double x, double (*fparams)[50], int n_strategy){
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
