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
#include <stdlib.h>

#define ARRAY_LEN 19

double nu;
double z;

form_factor_ctx ctx_spheres;
form_factor_ctx ctx_gdab;
form_factor_ctx ctx_broad_peak;

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
        387507.47872702,  371194.86504067,  349050.89244560,  
        322866.82508137,  293957.89067697,  263387.28834352,  
        232057.32254606,  200753.93208662,  170170.35793398,  
        140920.01392396,  113542.66309246,  88506.35995189,  
        66205.92725291,  46958.33549706,  30994.65355484,  
        18447.53465675,  9331.79386106,  3512.36330444,  643.29844390
    },
    .actual_spheres = { 0 }, 
    .expected_gdab = { 
        0.01314100,  0.00994096,  0.00744536,  0.00553396,  0.00408864,  
        0.00300617,  0.00220143,  0.00160669,  0.00116926,  0.00084882,  
        0.00061486,  0.00044454,  0.00032086,  0.00023125,  0.00016643,  
        0.00011964,  0.00008591,  0.00006162,  0.00004416
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
    double params_spheres[] = {10.0, 1.0};
    ctx_spheres.params = params_spheres;

    double params_gdab[] = {10.0, 0.5, 1e-4};
    ctx_gdab.params = params_gdab;

    double params_broad_peak[] = {10e5, 1000, 0.01, 2, 2};
    ctx_broad_peak.params = params_broad_peak;

    // setup the x (or r) array
    double r_array_spheres[ARRAY_LEN] = {
        1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0,  10.0,
        11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0
    };

    double r_array_gdab[ARRAY_LEN] = {
        15.,          18.54166667,  22.08333333,  25.625,       29.16666667,
        32.70833333,  36.25,        39.79166667,  43.33333333,  46.875,
        50.41666667,  53.95833333,  57.5,         61.04166667,  64.58333333,
        68.125,       71.66666667,  75.20833333,  78.75,             
    };

    double r_array_broad_peak[ARRAY_LEN] = {
        30.,   49.,   69.,   88.,   108.,  127.,  147.,  167.,  186.,  206.,  225., 
        225.,  245.,  265.,  284.,  304.,  323.,  343.,  362.,       
    };
    
    // COMPUTATIONS 
    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        z = r_array_spheres[i];
        hankel_transform_QWE_Key(
            nu, 
            form_factor_sphere, 
            z, 
            (void *)&ctx_spheres,
            &Gr[i], 
            250, 
            1e-9        
        );
        ctx.actual_spheres[i] = Gr[i]; 
    }
    
    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        z = r_array_gdab[i];
        hankel_transform_QWE_Key(
            nu, 
            form_factor_g_dab, 
            z, 
            (void *)&ctx_gdab,  
            &Gr[i],             
            250, 
            1e-9
        );
        ctx.actual_gdab[i] = Gr[i]; 
    }

    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        z = r_array_broad_peak[i];
        hankel_transform_QWE_Key(
            nu, 
            form_factor_broad_peak, 
            z, 
            (void *)&ctx_broad_peak, 
            &Gr[i],            
            150, 
            1e-9
        );
        ctx.actual_broad_peak[i] = Gr[i]; 
    }  
}


void tearDown(void) {}

void test_hankel_QWE_Key_regression_spheres(void) {
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

void test_hankel_QWE_Key_regression_gdab(void) {
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

void test_hankel_QWE_Key_regression_broad_peak(void) {
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

void test_hankel_QWE_Key_throws_error_when_not_converged(void) {
    /*
    Tests that hankel_QWE_Key throws expected
    error it doesn't converge
    */
    char captured[1024];
    int nu = 0;
    double z = 5.0;

    start_capture_stderr();
    int status = hankel_transform_QWE_Key(
        nu, 
        form_factor_g_dab, 
        z, 
        (void *)&ctx_gdab,
        &Gr[0], 
        4, 
        1e-9
    );
    stop_capture_stderr(captured, sizeof(captured));
    TEST_ASSERT_EQUAL_INT_MESSAGE(-4, status, "");
    TEST_ASSERT_EQUAL_STRING(
        captured, 
        "QWE_Key algorithm did not converge after maximum allowed intervals (4)\n"
    );
}

int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_hankel_QWE_Key_regression_spheres);
    RUN_TEST(test_hankel_QWE_Key_regression_gdab);
    RUN_TEST(test_hankel_QWE_Key_regression_broad_peak);
    RUN_TEST(test_hankel_QWE_Key_throws_error_when_not_converged);
    return UNITY_END();
}
