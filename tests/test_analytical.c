#include "src/utils/analytical_form_factors.h"

#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <stdbool.h>
#include "unity_config.h"
#include "unity.h"
#include <string.h>

#include "include/libhankel.h"
#include "src/utils/sasfit_integrate.h"


#define NUM_CASES 6
#define ARRAY_LEN 25


double nu = 0;
size_t int_strategy;
int z;
double params_spheres[50];
double params_gdab[50];

double r_array[ARRAY_LEN];
double G_analytic_spheres[ARRAY_LEN];
double G_analytic_gdab[ARRAY_LEN];

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
        -1158.122324, -3754.355368, -7278.680313, -11446.001590, -16047.001694, 
        -20912.464277, -25898.783335, -30880.874068, -35748.399431, -40403.736286, 
        -44760.976969, -48745.628409, -52294.852434, -55358.201871, -57898.900761, 
        -59895.836725, -61346.651880, -62272.843054, -62729.468960, -62729.468960, 
        -62729.468960, -62729.468960, -62729.468960, -62729.468960, -62729.468960, 
    },  
    .expected_gdab = {
        -7.34656e-05, -0.000225217, -0.000418221, -0.000634643, -0.000863458, 
        -0.00109729, -0.00133104, -0.00156112, -0.00178503, -0.00200103, -0.00220798, 
        -0.00240515, -0.00259213, -0.00276877, -0.00293509, -0.00309126, -0.00323752, 
        -0.00337421, -0.00350172, -0.00362046, -0.00373087, -0.00383338, -0.00392845, 
        -0.00401652, -0.00409801, 
    },  
};


void setUp(void) {
    /*
    Function that runs before tests, for set-up.    
    */
    
    // sphere params
    params_spheres[0] = 10; 
    params_spheres[1] = 1;

    params_gdab[0] = 10.0; 
    params_gdab[1] = 0.5;
    params_gdab[2] = 1e-4;

    // setup the x (or r) array
    double r_array[ARRAY_LEN] = {
        1.0,  2.0,  3.0,  4.0,  5.0,
        6.0,  7.0,  8.0,  9.0, 10.0,
        11.0, 12.0, 13.0, 14.0, 15.0,
        16.0, 17.0, 18.0, 19.0, 20.0,
        21.0, 22.0, 23.0, 24.0, 25.0
    };

    compute_analytical_spheres(&params_spheres, r_array, G_analytic_spheres, ARRAY_LEN);
    memcpy(ctx.actual_spheres, G_analytic_spheres, sizeof ctx.actual_spheres);

    compute_analytical_gdab(&params_gdab, r_array, G_analytic_gdab, ARRAY_LEN);
    memcpy(ctx.actual_gdab, G_analytic_gdab, sizeof ctx.actual_gdab);

    // printf("Analitycal solution:   ");
    // for (size_t i = 0; i < ARRAY_LEN; i++) {
    //     printf("%.6g, ", (G_analytic_gdab[i]));
    // };
    // printf("\n");
}


void tearDown(void) {}

void test_analytical_spheres(void) {
    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        TEST_ASSERT_DOUBLE_WITHIN(1e-4, ctx.expected_spheres[i], ctx.actual_spheres[i]);
    }
}

void test_analytical_gdab(void) {
    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        TEST_ASSERT_DOUBLE_WITHIN(1e-4, ctx.expected_gdab[i], ctx.actual_gdab[i]);
    }
}

int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_analytical_spheres);
    RUN_TEST(test_analytical_gdab);
    return UNITY_END();
}
