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
        387507.478876507, 371194.865559785, 349050.878110533, 322866.912359191, 
        293957.890774645, 263387.284336425, 232057.309283767, 200753.920990995, 
        170170.357980906, 140920.013857726, 113542.663215219, 88506.3597290707, 
        66205.9256771091, 46958.3386835007, 30994.6269460026, 18447.4103664975, 
        9331.51981974001, 3512.03954686726, 643.54697816898, -0.419462003739866, 
        0.671565816397894, 0.408622932460936, -0.572610815496937, 0.444648150293901, 
        -0.730794091442533, 

    },
    .actual_spheres = { 0 }, 
    .expected_gdab = { 
        0.0131409968968235, 0.0099409587529311, 0.00744536154206769, 0.00553395596369854, 
        0.00408864031395487, 0.00300616660224412, 0.00220142925163886, 0.00160668986388701,
        0.00116925997763928, 0.000848815153527063, 0.000614860095841608, 0.000444542888287165, 
        0.00032086242301512, 0.000231243932588246, 0.000166431099098308, 0.00011963763470055, 
        8.59047980625098e-05, 6.16207867381733e-05, 4.41601904858169e-05, 3.16199665583075e-05, 
        2.26226721933264e-05, 1.61732894051951e-05, 1.15543192494369e-05, 8.24883138794311e-06, 
        5.88501976636725e-06,
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
    //printf("de ogata, Gr-G0  \n");
    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        z = r_array_spheres[i];
        hankel_transform_DE_Ooura(
            nu, form_factor_sphere, z, &params_spheres, &Gr[i], 250, 1e-9
        );
        //printf("%.15g, ", (Gr[i]));
        ctx.actual_spheres[i] = Gr[i]; 
    }

    
    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        z = r_array_gdab[i];
        hankel_transform_DE_Ooura(
            nu, form_factor_g_dab, z, &params_gdab, &Gr[i], 250, 1e-9
        );
        ctx.actual_gdab[i] = Gr[i]; 
    }

    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        z = r_array_broad_peak[i];
        hankel_transform_DE_Ooura(
            nu, form_factor_broad_peak, z, &params_broad_peak, &Gr[i], 250, 1e-9
        );
        ctx.actual_broad_peak[i] = Gr[i]; 
    }  
}


void tearDown(void) {}

void test_hankel_DE_Ooura_regression_spheres(void) {
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

void test_hankel_DE_Ooura_regression_gdab(void) {
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

void test_hankel_DE_Ooura_regression_broad_peak(void) {
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
    RUN_TEST(test_hankel_DE_Ooura_regression_spheres);
    RUN_TEST(test_hankel_DE_Ooura_regression_gdab);
    RUN_TEST(test_hankel_DE_Ooura_regression_broad_peak);
    return UNITY_END();
}
