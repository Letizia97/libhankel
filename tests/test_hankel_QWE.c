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

double params[50];
double r_array[ARRAY_LEN];

hankel_inputs inputs;
double Gr[ARRAY_LEN];

// For running the tests
typedef struct {
    double actual[ARRAY_LEN];
    double expected[ARRAY_LEN];
} TestContext;

TestContext ctx = {
    .expected = {
        387507.478727, 371194.865046, 349050.892446, 322866.825081, 293957.890695,  
        263387.288344, 232057.322544, 200753.921990, 170170.357893, 140920.013924,  
        113542.663056, 88506.359952, 66205.927253, 46958.335485, 30994.653555, 
        18447.535099, 9331.793855, 3512.362861, 643.298173, -1.000000, -1.000000, 
        -0.000000, -0.000000, -0.000000, -0.000000,  
    },
    .actual = { 0 }, 
};



void setUp(void) {
    /*
    Function that runs before tests, for a shared set-up.
    Anything computed inside it is available in tests
    because previously declared outside the function.        
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

    nu = 0;

    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        z = r_array[i];
        Gr[i] = hankel_transform_QWE_Chave(nu, form_factor_sphere, z, &params, 250, 1e-9);
        printf("%f,  ", Gr[i]);
        ctx.actual[i] = Gr[i]; 
    }
    printf("\n");

}


void tearDown(void) {}

void test_hankel_QWE_regression(void) {
    /*
    Regression test for DHT strategies.
    */
    int failures = 0;

    if (!arrays_close(ctx.actual, ctx.expected, ARRAY_LEN, 1e-4)) {
        failures++;
        UnityPrint("Arrays didn't match");
        UNITY_OUTPUT_CHAR('\n');
    }

    TEST_ASSERT_EQUAL_MESSAGE(
        0,
        failures,
        "One or more arrays were not close enough"
    );
}


int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_hankel_QWE_regression);
    return UNITY_END();
}
