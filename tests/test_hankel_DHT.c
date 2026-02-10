#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <stdbool.h>
#include "utils/unity_config.h"
#include "unity.h"
#include <string.h>

#include "include/libhankel.h"
#include "include/form_factors.h"
#include "src/utils/sasfit_integrate.h"
#include "utils/test_utils.h" 



#define NUM_CASES 6
#define ARRAY_LEN 25

double nu;
double z;
size_t int_strategy;

double params[50];
double r_array[ARRAY_LEN];

hankel_inputs inputs;
double Gr[ARRAY_LEN];

// For running the tests
typedef struct {
    double actual[NUM_CASES][ARRAY_LEN];
    double expected[NUM_CASES][ARRAY_LEN];
    double rel_tol[NUM_CASES];
} TestContext;

TestContext ctx = {
    .expected = {
        {387675.515633, 371029.163618, 348543.605156, 322628.820713, 293627.711493, 263468.462343, 232210.835522, 200754.953535, 169913.695764, 141328.476997, 113451.018576, 88153.668256, 66335.095338, 47307.971825, 30989.654717, 18104.319228, 9120.610116, 3724.839675, 1018.234929, -23.871215, -232.938244, -139.283971, -8.875854, 64.713377, 76.053645}, 
        {387531.433668, 371629.814830, 350577.880311, 323723.742098, 293175.971072, 264237.938007, 231235.385078, 201256.966994, 170760.496514, 140310.011712, 113170.459667, 89103.269186, 66952.136369, 47029.915949, 30363.746456,  17880.491866, 9299.326596, 4102.411688, 1313.413421, 101.856515, -298.590691, -296.338187, -156.670157, -47.658432, 85.818980}, 
        {388886.903270, 371441.730622, 348775.029970, 322752.509972, 293781.898095, 263288.320330, 232351.232213, 200643.273981, 169986.161797, 141341.529033, 113382.932128, 88162.582887, 66380.832756, 47314.544723, 30963.287957, 18085.457896, 9123.270818, 3738.632396, 1028.777218, -22.863803, -238.232169, -145.489845, -12.293362, 65.199294, 78.939515}, 
        {386231.572315, 371450.709954, 348902.953186, 322889.532672, 293743.641156, 263296.956789, 232052.850582, 200942.032532, 169923.537891, 141095.964428, 113591.717101, 88262.287683, 66340.512223, 47111.388336, 30889.782673, 18292.841614, 9335.058230, 3693.883247, 796.476564, -193.268303, -207.517256, 30.243185, 156.011672, 117.369648, 3.325143}, 
        {387535.816919, 371150.250733, 349083.303256, 322942.714337, 293855.488835, 263362.120166, 232015.344905, 200669.147992, 170237.066265, 140885.343769, 113577.719928, 88445.145227, 66293.480530, 46885.601041, 30993.085203, 18526.810152, 9270.558881, 3462.287166, 747.950973, -30.619668, -55.379685, 20.376510, 27.348326, 0.526402, -13.846825}, 
        {387470.776802, 371207.185676, 349004.496127, 322909.376091, 293952.908042, 263402.816784, 232029.891423, 200730.644505, 170192.723764, 140960.942655, 113483.273545, 88560.754208, 66167.223872, 46970.971914, 31002.722413, 18453.786743, 9297.549934, 3498.030048, 738.887994, -48.654148, -60.971902, 25.124492, 34.597789, 3.200237, -16.795353}
    },
    .actual = { { 0 } },  // initializes expected[0][0], zeros everything else
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
    int_strategy = 6;
    for (size_t j = 0; j < NUM_CASES; ++j) {

        for (size_t i = 0; i < ARRAY_LEN; ++i) {
            z = r_array[i];
            Gr[i] = hankel_transform_DHT(nu, form_factor_sphere, z, &params, int_strategy);
            printf("%f,  ", Gr[i]);
            ctx.actual[j][i] = Gr[i]; 
        }
        printf("\n");
        int_strategy = int_strategy+1;
    }
}


void tearDown(void) {}

void test_hankel_DHT_regression(void) {
    /*
    Regression test for DHT strategies.
    */
    int failures = 0;

    for (size_t i = 0; i < NUM_CASES; ++i) {
        if (!arrays_close(ctx.actual[i], ctx.expected[i], ARRAY_LEN, 1e-4)) {
            failures++;

            UnityPrint("Array case failed: ");
            UnityPrintNumberUnsigned(i);
            UNITY_OUTPUT_CHAR('\n');
        }
    }

    TEST_ASSERT_EQUAL_MESSAGE(
        0,
        failures,
        "One or more arrays were not close enough"
    );
}

void test_hankel_DHT_throws_error_when_nu_equal_2(void) {
    /*
    Tests that hankel_transform_DHT throws expected
    error when nu is not zero nor 1.
    */
    nu = 2;
    z = 10.0;
    int_strategy = 7;
    char captured[1024];

    start_capture_stderr();
    int res = hankel_transform_DHT(nu, form_factor_sphere, z, &params, int_strategy);
    stop_capture_stderr(captured, sizeof(captured));
    TEST_ASSERT_EQUAL_INT_MESSAGE(-1, res, "");
    TEST_ASSERT_EQUAL_STRING(captured, "nu needs to be 0 or 1 in order to use the selected strategy\n");
}

void test_hankel_DHT_throws_error_when_int_strategy_wrong(void) {
    /*
    Tests that hankel_transform_DHT throws expected
    error when when int_strategy is not between 6 and 11.
    */
    nu = 0;
    z = 10.0;
    int_strategy = 1;
    char captured[1024];

    start_capture_stderr();
    int res = hankel_transform_DHT(nu, form_factor_sphere, z, &params, 1);
    stop_capture_stderr(captured, sizeof(captured));
    TEST_ASSERT_EQUAL_INT_MESSAGE(-1, res, "");
    TEST_ASSERT_EQUAL_STRING(captured, "Strategy number must be integer between 6 and 11\n");
}

int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_hankel_DHT_regression);
    RUN_TEST(test_hankel_DHT_throws_error_when_nu_equal_2);
    RUN_TEST(test_hankel_DHT_throws_error_when_int_strategy_wrong);
    return UNITY_END();
}
