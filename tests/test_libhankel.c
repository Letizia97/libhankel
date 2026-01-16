#include "tests/analytical_form_factors.h"

#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include "include/libhankel.h"
#include "src/utils/sasfit_integrate.h"


// gcc -Iinclude test/test_libhankel.c libhankel.a -o test_libhankel
// gcc -Iinclude test/test_libhankel.c libhankel.a -lgsl -lgslcblas -lm -o test_libhankel

int main() {
    double nu = 0;
    // double q = 0.1;
    double params[50] = {0};
    // int int_strategy = 6; // 6; //1;

    // g_dab params
    // params[0] = 10.0;
    // params[1] = 0.5;   
    // params[2] = 1e-4;

    // sphere params
    params[0] = 10; 
    params[1] = 1;

    double r_array[25] = {
        1.0,  2.0,  3.0,  4.0,  5.0,
        6.0,  7.0,  8.0,  9.0, 10.0,
        11.0, 12.0, 13.0, 14.0, 15.0,
        16.0, 17.0, 18.0, 19.0, 20.0,
        21.0, 22.0, 23.0, 24.0, 25.0
    };
    double G_analytic[25];
    double Gr[25];
    
    // Compute analytical
    compute_analytical_spheres(&params, r_array, G_analytic, 25);

    double ff_at_0(double q, hankel_inputs *params){
        return q * form_factor_sphere(q, params->fparams);
    }

    hankel_inputs inputs;
    inputs.function = form_factor_sphere;
	inputs.other_inputs[0] = nu;
	inputs.other_inputs[1] = 0;
    inputs.fparams=params;

    // Compute G0
    double G0 = sasfit_integrate_ctm(
        0,
        INFINITY,
        ff_at_0,
        &inputs,
        10000,
        0.001,
        1e-10
    );

    // Compute Hankel at r
    for (int i = 0; i < 25; i++) {
        int z = i+1;
        Gr[i] = hankel_transform_DE_Quadrature(nu, form_factor_sphere, z, &params, 250, 1e-9);
        printf("%f, ", (Gr[i] - G0) / (2 * M_PI));
        // hankel_transform_DHT(nu, form_factor_sphere, i, &params, int_strategy));
        // hankel_transform_FBT(nu, form_factor_sphere, i, &params, int_strategy, 250, 1e-3));
        // hankel_transform_DE_Ogata(nu, form_factor_sphere, i, &params, 250, 1e-3));
        // hankel_transform_QWE_Key(nu, form_factor_sphere, i, &params, 350, 1e-9));
    };
    printf("\n");

    for (int i = 0; i < 25; i++) {
        printf("%f, ", (G_analytic[i]));
    };
    
    return 0;
}
