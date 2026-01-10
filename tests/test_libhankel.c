
#include <stdio.h>
#include "include/libhankel.h"


// gcc -Iinclude test/test_libhankel.c libhankel.a -o test_libhankel
// gcc -Iinclude test/test_libhankel.c libhankel.a -lgsl -lgslcblas -lm -o test_libhankel

int main() {
    double nu = 0;
    double q = 0.1;
    double params[50] = {0};
    int int_strategy = 2; // 6; //1;

    // g_dab params
    // params[0] = 10.0;
    // params[1] = 0.5;   
    // params[2] = 1e-4;

    // sphere params
    params[0] = 10; 
    params[1] = 1;


    printf("Result of calling g_dab: %f \n", form_factor_sphere(q, &params));
    
    printf("Result of calling hankel: \n");
    for (double i = 1; i < 25; i++) {
        // printf(" %f, ",  hankel_transform_no_params(nu, form_factor_sphere, i, &params, int_strategy));
        // printf(" %f, ",  hankel_transform_FBT(nu, form_factor_sphere, i, &params, int_strategy, 250, 1e-3));
        // printf(" %f, ",  hankel_transform_DE_Quadrature(nu, form_factor_sphere, i, &params, 250, 1e-9));
        // printf(" %f, ",  hankel_transform_DE_Ogata(nu, form_factor_sphere, i, &params, 250, 1e-3));
        printf(" %.17f, ",  hankel_transform_QWE_Key(nu, form_factor_sphere, i, &params, 350, 1e-9));
    };

    return 0;

}

// Seems that the fbt transform is not working correctly because i ge valeus different from those in python

