#include "src/utils/analytical_form_factors.h"

#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <stdbool.h>
#include "unity_config.h"
#include "unity.h"
#include <string.h>

#include "include/libhankel.h"
#include "include/form_factors.h"
#include "src/utils/sasfit_integrate.h"

#define NUM_CASES 6
#define ARRAY_LEN 25

double nu = 0;
size_t int_strategy;
int z;

size_t n_params_spheres = 2;
size_t n_params_gdab = 3;

double params_spheres[2];
double params_gdab[3];

double r_array_spheres[ARRAY_LEN];
double r_array_gdab[ARRAY_LEN];

double res_spheres[ARRAY_LEN];
double res_gdab[ARRAY_LEN];

hankel_inputs inputs;
double Gr[ARRAY_LEN];
double G0;
double G_sesans[ARRAY_LEN];

// For running the tests
typedef struct {
    double actual_spheres[ARRAY_LEN];
    double actual_gdab[ARRAY_LEN];
    double expected_spheres[ARRAY_LEN];
    double expected_gdab[ARRAY_LEN];
} TestContext;

TestContext ctx = {
    .actual_spheres = {0},
    .actual_gdab = {0},
    .expected_spheres = {
        9722.84094427162, 129.64612039143, 6.83094332765291, 28.9919993496471, 
        23.7834230065539, 10.9350274698768, 2.5466918177704, 0.0369989644552898, 
        0.504886496554671, 1.18806736999471, 1.07559381275325, 0.498837948141705, 
        0.0717102119680533, 0.0172437738000385, 0.154603937929802, 
        0.230000915953415, 0.165628310445537, 0.0530781871256424, 
        0.000451701803214946, 0.0238474803585097, 0.0637547626792619,
        0.0668300062717888, 0.0347749850414803, 0.00493014922872318,
        0.00242399844153717
    },  
    .expected_gdab = {
        0.0386509844824505, 0.033247733037704, 0.0059906836301716, 
        0.000413335562906358, 0.00421046656445929, 0.00301952811898612, 
        0.000270578715140846, 0.000372204882984811, 0.00107361158846175, 
        0.00051441098212875, 1.98106457281109e-06, 0.000240613262724929, 
        0.000358215907523952, 9.63920710994842e-05, 1.25705109248615e-05, 
        0.000145080403092356, 0.000128553677409393, 1.27141486164948e-05, 
        2.64407254593612e-05, 8.33460220045826e-05, 4.39973323194924e-05, 
        6.33002337339746e-08, 2.86891932932902e-05, 4.49989296028058e-05, 
        1.24491585407486e-05
    },  
};


void setUp(void) {
    /*
    Function that runs before tests, for set-up.    
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

    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        ctx.actual_spheres[i] = form_factor_sphere(
            r_array_spheres[i], 
            params_spheres,
            n_params_spheres
        );
        ctx.actual_gdab[i] = form_factor_sphere(
            r_array_gdab[i], 
            params_gdab,
            n_params_gdab
        );
        //printf("%.15g, ", ctx.actual_gdab[i]);
    }
}


void tearDown(void) {}

void test_form_factor_spheres(void) {
    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        TEST_ASSERT_DOUBLE_WITHIN(
            1e-4, 
            ctx.expected_spheres[i], 
            ctx.actual_spheres[i]
        );
    }
}

void test_form_factor_gdab(void) {
    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        TEST_ASSERT_DOUBLE_WITHIN(
            1e-4, 
            ctx.expected_gdab[i], 
            ctx.actual_gdab[i]
        );
    }
}

int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_form_factor_spheres);
    RUN_TEST(test_form_factor_gdab);
    return UNITY_END();
}
