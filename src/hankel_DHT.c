
#include "libhankel.h"

#include <math.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_bessel.h>

#include "../src/utils/strateg6_const.h"
#include "../src/utils/strateg7_const.h"
#include "../src/utils/strateg8_const.h"
#include "../src/utils/strateg9_const.h"
#include "../src/utils/strateg10_const.h"
#include "../src/utils/strateg11_const.h"

/*
This file contains functions corresponding to strategies 6-11 in SASfit 
(see inline comments for specific names).

They have been grouped together under one function, as they are very similar, and
changing the n_strategy parameter allows to switch between them
*/



/** 
 * @brief Computes Hankel transform, using digital filters.
 * @note Corresponds to strategies 6-11 in SASfit.
 * @note Does not allow any error control.
 * @note Perform wells for simple form factors, but it does struggle 
 *       with oscillatory ones. For oscillatory form factors, it is 
 *       advisable to start with this method for a rough Hankel transform
 *       computation, and then refine it with "hankel_transform_QWE_Chave"
 *       or "hankel_transform_QWE_Key".
 * 
 * @param nu         order of bessel function - must be 0 or 1
 * @param f          pointer to form factor function
 * @param x          value at which to compute the transform
 * @param f_params    params for form factor
 * @param output     pointer to var containing output from transform 
 * @param n_strategy an integer number between 6 and 11 
 *                   (determines which weightings to use for the transform
 *                   and corresponds to SASfit strategies 6-11)
 */
double hankel_transform_DHT(
    int nu, 
    double (*f)(double, double (*)[50]), 
    double x, 
    double (*f_params)[50], 
    double *output,
    int n_strategy) {

    double res = 0;
    double lambda;
    unsigned int i;
    unsigned int ind;

    ind = nu + 1;

    if (!(nu==0 || nu==1)) {
        fprintf(stderr, "nu needs to be 0 or 1 in order to use the selected strategy\n");
        return -1;
    }

    if (!(n_strategy >= 6 && n_strategy <= 11)) {
        fprintf(stderr, "Strategy number must be integer between 6 and 11\n");
        return -2;
    }

    switch(n_strategy){

        case 6:{
            // HANKEL_GUPTASARMA_97
            if (nu==0) {
                for (i=0; i<120; i++) {
                    lambda = pow(10.0E0, (aJ0 + i*sJ0)) / x;
                    res = res + (*f)(lambda, f_params) * lambda * WJ0[i] / x;
                }
            } else {
                for (i=0; i<140; i++) {
                    lambda = pow(10.0E0, (aJ1 + i*sJ1)) / x;
                    res = res + (*f)(lambda, f_params) * lambda * WJ1[i] / x;
                }
            }
            break;
        }
        case 7:{
            // HANKEL_GUPTASARMA_97_FAST
            if (nu==0) {
                for (i=0; i<61; i++) {
                    lambda = pow(10.0E0, (aJ0Fast + i*sJ0Fast))/x;
                    res = res+(*f)(lambda, f_params) * lambda * WJ0Fast[i] / x;
                }
            } else {
                for (i=0; i<47; i++) {
                    lambda = pow(10.0E0, (aJ1Fast + i*sJ1Fast)) / x;
                    res = res + (*f)(lambda, f_params) * lambda * WJ0Fast[i] / x;
                }
            }
            break;
        }
        case 8:{
            // HANKEL_KEY_51
            for (i=0; i<51; i++) {
                lambda = KK51Hankel[i][0] / x;
                res = res + (*f)(lambda, f_params) * lambda * KK51Hankel[i][ind] / x;
            }
            break;
        }
        case 9:{
            // HANKEL_KEY_101
            for (i=0; i<101; i++) {
                lambda = KK101Hankel[i][0]/x;
                res = res + (*f)(lambda, f_params) * lambda * KK101Hankel[i][ind] / x;
            }
            break;
        }
        case 10:{
            // HANKEL_KEY_201
            for (i=0; i<201; i++) {
                lambda = KK201Hankel[i][0] / x;
                res = res + (*f)(lambda, f_params) * lambda * KK201Hankel[i][ind] / x;
            }
            break;
        }
        case 11:{
            // HANKEL_ANDERSON_801
            for (i=0; i<801; i++) {
                lambda = WA801Hankel[i][0] / x;
                res = res + (*f)(lambda, f_params) * lambda * WA801Hankel[i][ind] / x;
            } 
            break;

        }

    }
    *output = res;
    return 0;
}
