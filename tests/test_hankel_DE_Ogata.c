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
        387485.814841457, 371194.098327246, 349051.627314929, 322866.701679252,
        293957.714319064, 263387.299263421, 232057.360059439, 200753.961948608, 
        170170.4100575, 140920.060308871, 113542.745812402, 88506.3737683421, 
        66206.0134182743, 46958.3822491802, 30994.6595042601, 18447.7120479117, 
        9332.22332872738, 3512.55955045874, 643.163700220762, -2.63907790416413, 
        -0.97139642061165, 0.900350878602035, 0.76008614871604, -0.577281648088827, 
        0.948387334421353, 
    },
    .actual_spheres = { 0 }, 
    .expected_gdab = { 
        0.0131410051769628, 0.00994097008786127, 0.00744537784848448, 
        0.00553397819706947, 0.00408866884933503, 0.00300620247657538, 
        0.00220147329769634, 0.00160674291367746, 0.00116932286096173, 
        0.000848888697775925, 0.000614945125807072, 0.000444640225954316, 
        0.000320972887362585, 0.000231368339386456, 0.000166570260707958, 
        0.000119792359871521, 8.60761786532446e-05, 6.18093625176832e-05, 
        4.43668947863511e-05, 3.18454671319162e-05, 2.2867758083941e-05, 
        1.64388080641872e-05, 1.18409888712984e-05, 8.5574251807065e-06, 
        6.21630545540752e-06, 
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

    // Read broad peak expected values from file
    const char *filename = "../tests/data/broadpeakHT.txt";
    const char *column   = "xi1000";
    size_t n             = sizeof(r_array_broad_peak) / sizeof(r_array_broad_peak[0]);
    read_values_by_rows(
        filename, 
        column, 
        r_array_broad_peak, 
        n, 
        ctx.expected_broad_peak
    );

    
    // COMPUTATIONS 
    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        z = r_array_spheres[i];
        hankel_transform_DE_Ogata(
            nu, form_factor_sphere, z, &params_spheres, &Gr[i], 150, 1e-3
        );

        ctx.actual_spheres[i] = Gr[i]; 
    }

    //printf("de ogata, Gr-G0  \n");
    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        z = r_array_gdab[i];
        hankel_transform_DE_Ogata(
            nu, form_factor_g_dab, z, &params_gdab, &Gr[i], 150, 1e-3
        );
        //printf("%.15g, ", (Gr[i]));
        ctx.actual_gdab[i] = Gr[i]; 
    }

    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        z = r_array_broad_peak[i];
        hankel_transform_DE_Ogata(
            nu, form_factor_broad_peak, z, &params_broad_peak, &Gr[i], 150, 0.2e-3
        );
        ctx.actual_broad_peak[i] = Gr[i]; 
    }  
}


void tearDown(void) {}

void test_hankel_DE_Ogata_regression_spheres(void) {
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

void test_hankel_DE_Ogata_regression_gdab(void) {
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

void test_hankel_DE_Ogata_regression_broad_peak(void) {
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

int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_hankel_DE_Ogata_regression_spheres);
    RUN_TEST(test_hankel_DE_Ogata_regression_gdab);
    // probably need to replace this with an actual regression
    //RUN_TEST(test_hankel_DE_Ogata_regression_broad_peak);
    return UNITY_END();
}
