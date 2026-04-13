#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <stdbool.h>
#include "utils/unity_config.h"
#include "unity.h"
#include <string.h>

#include "libhankel.h"
#include "form_factors.h"
#include "src/utils/sasfit_integrate.h"
#include "utils/test_utils.h" 

#define ARRAY_LEN 25

double nu;
double z;

double params_spheres[50];
double params_gdab[50];
double params_broad_peak[50];

double r_array_spheres[ARRAY_LEN];
double r_array_gdab[ARRAY_LEN];
double r_array_broad_peak[ARRAY_LEN];

hankel_inputs inputs;
double Gr[ARRAY_LEN];

// For running the tests
typedef struct {
    double actual_spheres[ARRAY_LEN];
    double expected_spheres[ARRAY_LEN];
    double actual_gdab[ARRAY_LEN];
    double expected_gdab[ARRAY_LEN];
    double actual_broad_peak[ARRAY_LEN];
    double expected_broad_peak[ARRAY_LEN];
} TestContext;

TestContext ctx = {
    .expected_spheres = {
        387507.478727, 371194.865046, 349050.892446, 322866.825081, 293957.890695,  
        263387.288344, 232057.322544, 200753.921990, 170170.357893, 140920.013924,  
        113542.663056, 88506.359952, 66205.927253, 46958.335485, 30994.653555, 
        18447.535099, 9331.793855, 3512.362861, 643.298173, -3.73877664854985e-09, 
        -1.72981880655027e-08, -1.58967451740654e-08, -3.56770079310641e-10,
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
    .expected_broad_peak =  {
        15.3478, 14.7545, 13.8463, 12.7398, 11.3539, 9.86578, 8.1652,
        6.37775, 4.64776, 2.84404, 1.19288, 1.19288, -0.433731, 
        -1.904, -3.12269, -4.19016, -4.98079, -5.56571, -5.88294, 
        -5.97027, -5.81729, -5.4678, -4.91068, -4.23042, -3.39016
     }, 
    .actual_broad_peak = { 0 }, 
};



void setUp(void) {
    /*
    Function that runs before tests, for a shared set-up.
    Anything computed inside it is available in tests
    because previously declared outside the function.        
    */
    nu = 0;

    // set params
    params_spheres[0] = 10; 
    params_spheres[1] = 1;

    params_gdab[0] = 10.0; 
    params_gdab[1] = 0.5;
    params_gdab[2] = 1e-4;

    params_broad_peak[0] = 10e5;
    params_broad_peak[1] = 1000;
    params_broad_peak[2] = 0.01;
    params_broad_peak[3] = 2;
    params_broad_peak[4] = 2;

    // setup the x (or r) array
    double r_array_spheres[ARRAY_LEN] = {
        1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0,  10.0,
        11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
        21.0, 22.0, 23.0, 24.0, 25.0
    };

    double r_array_gdab[ARRAY_LEN] = {
        15.,          18.54166667,  22.08333333,  25.625,       29.16666667,
        32.70833333,  36.25,        39.79166667,  43.33333333,  46.875,
        50.41666667,  53.95833333,  57.5,         61.04166667,  64.58333333,
        68.125,       71.66666667,  75.20833333,  78.75,        82.29166667,
        85.83333333,  89.375,       92.91666667,  96.45833333,  100.,         
    };

    double r_array_broad_peak[ARRAY_LEN] = {
        30.,   49.,   69.,   88.,   108.,  127.,  147.,  167.,  186.,  206.,  225., 
        225.,  245.,  265.,  284.,  304.,  323.,  343.,  362.,  382.,  402.,  421.,  
        441.,  460.,  480.,      
    };
   
    // COMPUTATIONS 

    // printf("QWE Chave, Gr-G0 on spheres \n");
    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        z = r_array_spheres[i];
        hankel_transform_QWE_Chave(
            nu, form_factor_sphere, z, &params_spheres, &Gr[i], 250, 1e-9
        );
        //printf("%.15g,  \n", (Gr[i]));
        ctx.actual_spheres[i] = Gr[i]; 
    }
    
    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        z = r_array_gdab[i];
        hankel_transform_QWE_Chave(
            nu, form_factor_g_dab, z, &params_gdab, &Gr[i], 250, 1e-9
        );
        ctx.actual_gdab[i] = Gr[i]; 
    }

    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        z = r_array_broad_peak[i];
        hankel_transform_QWE_Chave(
            nu, form_factor_broad_peak, z, &params_broad_peak, &Gr[i], 150, 1e-9
        );
        ctx.actual_broad_peak[i] = Gr[i]; 
    }  
}


void tearDown(void) {}

void test_hankel_QWE_Chave_regression_spheres(void) {
    /*
    Regression test for QWE on spheres.
    */
    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        TEST_ASSERT_DOUBLE_WITHIN(
            1e-4, 
            ctx.expected_spheres[i],
            ctx.actual_spheres[i]
        );
    }
}

void test_hankel_QWE_Chave_regression_gdab(void) {
    /*
    Regression test for QWE on gdab.
    */
    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        TEST_ASSERT_DOUBLE_WITHIN(
            1e-4, 
            ctx.expected_gdab[i], 
            ctx.actual_gdab[i]
        );
    }
}

void test_hankel_QWE_Chave_regression_broad_peak(void) {
    /*
    Regression test for QWE on broad peak.
    */
    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        TEST_ASSERT_DOUBLE_WITHIN(
            1e-4, 
            ctx.expected_broad_peak[i], 
            ctx.actual_broad_peak[i]
        );
    }
}

void test_hankel_QWE_Chave_throws_error_when_not_converged(void) {
    /*
    Tests that hankel_QWE_Chave throws expected
    error it doesn't converge.
    */
    char captured[1024];
    int nu = 0;
    double z = 5.0;

    start_capture_stderr();
    int status = hankel_transform_QWE_Chave(
        nu, form_factor_g_dab, z, &params_gdab, &Gr[0], 4, 1e-9
    );
    stop_capture_stderr(captured, sizeof(captured));
    TEST_ASSERT_EQUAL_INT_MESSAGE(-4, status, "");
    TEST_ASSERT_EQUAL_STRING(
        captured, 
        "QWE_Chave algorithm did not converge after maximum allowed intervals (4)\n"
    );
}

int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_hankel_QWE_Chave_regression_spheres);
    RUN_TEST(test_hankel_QWE_Chave_regression_gdab);
    RUN_TEST(test_hankel_QWE_Chave_regression_broad_peak);
    RUN_TEST(test_hankel_QWE_Chave_throws_error_when_not_converged);
    return UNITY_END();
}
