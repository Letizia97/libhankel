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
        387507.479804014, 371194.889604113, 349050.937049571,  
        322866.922336953, 293957.907818019, 263387.963887526,  
        232057.193789323, 200754.164640484, 170171.104078434,  
        140922.587792766, 113545.542472399, 88505.3243843713,  
        66204.9013454499, 46966.5166892448, 30987.387143271,  
        18459.4870999164, 9333.31546877299, 3491.19343376075,  
        699.480219589015,   
    },
    .actual_spheres = { 0 }, 
    .expected_gdab = { 
        0.013141156724346, 0.0099412012226085, 0.00744570498637932,  
        0.00553441753202811, 0.00408923631685086, 0.00300691371816697,  
        0.00220234362649212, 0.00160778728125935, 0.00117055582558994,  
        0.000850324394106815, 0.000616597236201586, 0.000446521953475126,  
        0.00032309693054483, 0.000233746868850275, 0.000169214897682618,  
        0.000122714156685465, 8.92856012262321e-05, 6.53162748940863e-05,  
        4.81805457319087e-05,  
    },
    .actual_gdab = { 0 }, 
    .expected_broad_peak =  {
        15.3241805903474, 14.7749665280258, 13.8451152654997,  
        12.7390232840474, 11.3539775134826, 9.86576027958892,  
        8.1651795673557, 6.37775908103334, 4.64776019501934,  
        2.84404211244366, 1.1928794134434, 1.1928794134434,  
        -0.433730756150041, -1.90399578725017, -3.12268829208393,  
        -4.19016310517701, -4.98078692314085, -5.56571377812195,  
        -5.88293888026089,  
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
        hankel_transform_FBT(
            nu, form_factor_sphere, z, &params_spheres, &Gr[i], 0, 250, 1e-3
        );
        ctx.actual_spheres[i] = Gr[i]; 
    }
    
    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        z = r_array_gdab[i];
        hankel_transform_FBT(
            nu, form_factor_g_dab, z, &params_gdab, &Gr[i], 0, 250, 1e-3
        );
        ctx.actual_gdab[i] = Gr[i]; 
    }

    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        z = r_array_broad_peak[i];
        hankel_transform_FBT(
            nu, form_factor_broad_peak, z, &params_broad_peak, &Gr[i], 0, 250, 1e-3
        );
        printf("%.15g,  \n", (Gr[i]));
        ctx.actual_broad_peak[i] = Gr[i]; 
    }  
}


void tearDown(void) {}

void test_hankel_FBT_regression_spheres(void) {
    /*
    Regression test for FBT on spheres.
    */
    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        TEST_ASSERT_DOUBLE_WITHIN(
            1e-4, 
            ctx.expected_spheres[i], 
            ctx.actual_spheres[i]
        );
    }
}

void test_hankel_FBT_regression_gdab(void) {
    /*
    Regression test for FBT on gdab.
    */
    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        TEST_ASSERT_DOUBLE_WITHIN(
            1e-4, 
            ctx.expected_gdab[i], 
            ctx.actual_gdab[i]
        );
    }
}

void test_hankel_FBT_regression_broad_peak(void) {
    /*
    Regression test for FBT on broad peak.
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
    RUN_TEST(test_hankel_FBT_regression_spheres);
    RUN_TEST(test_hankel_FBT_regression_gdab);
    RUN_TEST(test_hankel_FBT_regression_broad_peak);
    return UNITY_END();
}
