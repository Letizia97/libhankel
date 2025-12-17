
#include <stdio.h>
#include "include/libhankel.h"


// gcc -Iinclude test/test_libhankel.c libhankel.a -o test_libhankel
//  gcc -Iinclude test/test_libhankel.c libhankel.a -lgsl -lgslcblas -lm -o test_libhankel

int main() {
    double nu = 0;
    double q = 0.1;
    double params[50] = {0};
    int int_strategy = 1;//6; //1;

    // g_dab params
    // params[0] = 10.0;
    // params[1] = 0.5;   
    // params[2] = 1e-4;

    // sphere params
    params[0] = 100; //1e-9;
    params[1] = 1e-4;

    // printf("Result of calling g_dab: %f \n", form_factor_sphere(q, &params));

    printf("Result of calling g_dab: %f \n", form_factor_sphere(q, &params));


    printf("Result of calling hankel: \n");
    for (double i = 10; i < 100; i++) {
        //printf(" %f, ",  hankel_transform_no_params(nu, form_factor_sphere, i, &params, int_strategy));
        printf(" %f, ",  compute_hankel_FBT(nu, form_factor_sphere, i, &params, int_strategy, 150, 1e-3));
    };

    return 0;

}

