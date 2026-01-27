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
double params[50];
double r_array[ARRAY_LEN];
double G_analytic[ARRAY_LEN];

hankel_inputs inputs;
double Gr[ARRAY_LEN];
double G0;
double G_sesans[ARRAY_LEN];

// For running the tests
typedef struct {
    double actual[ARRAY_LEN];
    double expected[ARRAY_LEN];
} TestContext;

TestContext ctx = {
    .actual = {0},
    .expected = {
        -1158.122324, -3754.355368, -7278.680313, -11446.001590, -16047.001694, 
        -20912.464277, -25898.783335, -30880.874068, -35748.399431, -40403.736286, 
        -44760.976969, -48745.628409, -52294.852434, -55358.201871, -57898.900761, 
        -59895.836725, -61346.651880, -62272.843054, -62729.468960, -62729.468960, 
        -62729.468960, -62729.468960, -62729.468960, -62729.468960, -62729.468960, 
    },  
};


void setUp(void) {
    /*
    Function that runs before tests, for set-up.    
    */
    
    // sphere params
    params[0] = 10; 
    params[1] = 1;

    // setup the x (or r) array
    double r_array[ARRAY_LEN] = {
        1.0,  2.0,  3.0,  4.0,  5.0,
        6.0,  7.0,  8.0,  9.0, 10.0,
        11.0, 12.0, 13.0, 14.0, 15.0,
        16.0, 17.0, 18.0, 19.0, 20.0,
        21.0, 22.0, 23.0, 24.0, 25.0
    };
    
    // Compute analytical
    compute_analytical_spheres(&params, r_array, G_analytic, ARRAY_LEN);
    memcpy(ctx.actual, G_analytic, sizeof ctx.actual);

    // printf("Analitycal solution:   ");
    // for (size_t i = 0; i < ARRAY_LEN; i++) {
    //     printf("%f, ", (G_analytic[i]));
    // };
    // printf("\n");
}


void tearDown(void) {}

void test_analytical(void) {
    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        TEST_ASSERT_DOUBLE_WITHIN(1e-4, ctx.expected[i], ctx.actual[i]);
    }
}

int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_analytical);
    return UNITY_END();
}
