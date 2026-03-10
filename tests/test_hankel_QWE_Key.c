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

#define ARRAY_LEN 19

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
    .expected_broad_peak =  { 0 }, 
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

    // Read broad peak expected values from file
    const char *filename = "../tests/data/broadpeakHT.txt";
    const char *column   = "xi1000";
    size_t n             = sizeof(r_array_broad_peak) / sizeof(r_array_broad_peak[0]);
    read_values_by_rows(filename, column, r_array_broad_peak, n, ctx.expected_broad_peak);

    
    // COMPUTATIONS 
    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        z = r_array_spheres[i];
        hankel_transform_QWE_Key(
            nu, form_factor_sphere, z, &params_spheres, &Gr[i], 250, 1e-9
        );
        ctx.actual_spheres[i] = Gr[i]; 
    }
    
    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        z = r_array_gdab[i];
        hankel_transform_QWE_Key(
            nu, form_factor_g_dab, z, &params_gdab, &Gr[i], 250, 1e-9
        );
        ctx.actual_gdab[i] = Gr[i]; 
    }

    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        z = r_array_broad_peak[i];
        hankel_transform_QWE_Key(
            nu, form_factor_broad_peak, z, &params_broad_peak, &Gr[i], 150, 1e-9
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
        TEST_ASSERT_DOUBLE_WITHIN(1e-4, ctx.expected_spheres[i], ctx.actual_spheres[i]);
    }
}

void test_hankel_QWE_Key_regression_gdab(void) {
    /*
    Regression test for QWE on gdab.
    */
    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        TEST_ASSERT_DOUBLE_WITHIN(1e-4, ctx.expected_gdab[i], ctx.actual_gdab[i]);
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

void test_hankel_QWE_Key_throws_error_when_nu_wrong(void) {
    /*
    Tests that hankel_QWE_Key throws expected
    error when nu is neither 0 nor 1.
    */
    char captured[1024];
    int nu = 2;
    double z = 5.0;

    start_capture_stderr();
    int status = hankel_transform_QWE_Key(
        nu, form_factor_g_dab, z, &params_gdab, &Gr[0], 250, 1e-9
    );
    stop_capture_stderr(captured, sizeof(captured));
    TEST_ASSERT_EQUAL_INT_MESSAGE(-1, status, "");
    TEST_ASSERT_EQUAL_STRING(
        captured, 
        "nu needs to be 0 or 1 in order to use QWE_Key\n"
    );
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
        nu, form_factor_g_dab, z, &params_gdab, &Gr[0], 4, 1e-9
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
    RUN_TEST(test_hankel_QWE_Key_throws_error_when_nu_wrong);
    RUN_TEST(test_hankel_QWE_Key_throws_error_when_not_converged);
    return UNITY_END();
}
