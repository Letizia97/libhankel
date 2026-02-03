#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <stdbool.h>
#include "utils/unity_config.h"
#include "unity.h"
#include <string.h>

#include "include/libhankel.h"
#include "src/utils/sasfit_integrate.h"
#include "utils/test_utils.h" 

#define ARRAY_LEN 25
#define MAX_COLS 25
#define MAX_ROWS 6

double nu;
double z;

double params_spheres[50];
double params_gdab[50];

double r_array_spheres[ARRAY_LEN];
double r_array_gdab[ARRAY_LEN];

hankel_inputs inputs;
double Gr[ARRAY_LEN];
double G0_spheres;

// For running the tests
typedef struct {
    double actual_spheres[ARRAY_LEN];
    double expected_spheres[ARRAY_LEN];
    double actual_gdab[ARRAY_LEN];
    double expected_gdab[ARRAY_LEN];
} TestContext;

TestContext ctx = {
    .expected_spheres = {
        387507.478727, 371194.865046, 349050.892446, 322866.825081, 293957.890695,  
        263387.288344, 232057.322544, 200753.921990, 170170.357893, 140920.013924,  
        113542.663056, 88506.359952, 66205.927253, 46958.335485, 30994.653555, 
        18447.535099, 9331.793855, 3512.362861, 643.298173, -1.000000, -1.000000, 
        -0.000000, -0.000000, -0.000000, -0.000000,  
    },
    .actual_spheres = { 0 }, 
    .expected_gdab = { 
        0.013140998,  0.0099409587,  0.0074453617,  0.0055339564,  0.0040886406,  
        0.003006167,  0.0022014298,  0.0016066905,  0.0011692607,  0.00084881602,  
        0.00061486109,  0.00044454403,  0.00032086372,  0.0002312454,  0.00016643275,  
        0.00011963947,  8.5907116e-05,  6.1623338e-05,  4.4163123e-05,  3.1623165e-05, 
        2.2626149e-05,  1.617712e-05,  1.1558453e-05,  8.2532795e-06,  5.8897925e-06,  
    },
    .actual_gdab = { 0 }, 
};

// Needed for computing G0
// double spheres_ff_at_0(double q, hankel_inputs *params){
//     return q * form_factor_sphere(q, params->fparams);
// }


void setUp(void) {
    /*
    Function that runs before tests, for a shared set-up.
    Anything computed inside it is available in tests
    because previously declared outside the function.        
    */
    
    // set params
    params_spheres[0] = 10; 
    params_spheres[1] = 1;

    params_gdab[0] = 10.0; 
    params_gdab[1] = 0.5;
    params_gdab[2] = 1e-4;

    // setup the x (or r) array
    double r_array_spheres[ARRAY_LEN] = {
        1.0,  2.0,  3.0,  4.0,  5.0,
        6.0,  7.0,  8.0,  9.0, 10.0,
        11.0, 12.0, 13.0, 14.0, 15.0,
        16.0, 17.0, 18.0, 19.0, 20.0,
        21.0, 22.0, 23.0, 24.0, 25.0
    };

    double r_array_gdab[ARRAY_LEN] = {
        15.,          18.54166667,  22.08333333,  25.625,       29.16666667,
        32.70833333,  36.25,        39.79166667,  43.33333333,  46.875,
        50.41666667,  53.95833333,  57.5,         61.04166667,  64.58333333,
        68.125,       71.66666667,  75.20833333,  78.75,        82.29166667,
        85.83333333,  89.375,       92.91666667,  96.45833333,  100.,         
    };


    
    nu = 0;

    // Compute G0
    // inputs.function = form_factor_sphere;
	// inputs.other_inputs[0] = nu;
	// inputs.other_inputs[1] = 0;
    // inputs.fparams=params_spheres;

    // G0_spheres = sasfit_integrate_ctm(
    //     0,
    //     INFINITY,
    //     spheres_ff_at_0,
    //     &inputs,
    //     10000,
    //     0.001,
    //     1e-20
    // );
    // printf("&&&&&&&&&&&&&&&&&     G0 %f, ", G0_spheres);



    // printf("QWE Chave, Gr-G0 on spheres \n");
    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        z = r_array_spheres[i];
        Gr[i] = hankel_transform_QWE_Chave(nu, form_factor_sphere, z, &params_spheres, 250, 1e-9);
        // printf("%.15g,  \n", (Gr[i]-G0_spheres)/ (2 * M_PI));
        ctx.actual_spheres[i] = Gr[i]; 
    }
    // printf("\n");

    
    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        z = r_array_gdab[i];
        Gr[i] = hankel_transform_QWE_Chave(nu, form_factor_g_dab, z, &params_gdab, 250, 1e-9);
        ctx.actual_gdab[i] = Gr[i]; 
    }

}


void tearDown(void) {}

void test_hankel_QWE_regression_spheres(void) {
    /*
    Regression test for DHT strategies.
    */
    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        TEST_ASSERT_DOUBLE_WITHIN(1e-4, ctx.expected_spheres[i], ctx.actual_spheres[i]);
    }
}

void test_hankel_QWE_regression_gdab(void) {
    /*
    Regression test for DHT strategies.
    */
    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        TEST_ASSERT_DOUBLE_WITHIN(1e-4, ctx.expected_gdab[i], ctx.actual_gdab[i]);
    }
}


int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_hankel_QWE_regression_spheres);
    RUN_TEST(test_hankel_QWE_regression_gdab);
    return UNITY_END();
}
