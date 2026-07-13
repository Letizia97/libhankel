#include "src/utils/analytical_form_factors.h"

// clang-format off
#include "unity_config.h"
#include "unity.h"
// clang-format on
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "include/libhankel.h"
#include "src/utils/sasfit_integrate.h"

#define NUM_CASES 6
#define ARRAY_LEN 25

int nu = 0;
int int_strategy;
int z;
double params_spheres[50];
double params_gdab[50];

double r_array_spheres[ARRAY_LEN];
double r_array_gdab[ARRAY_LEN];

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
    .expected_spheres =
        {
            -1158.12232352,  -3754.35536762,  -7278.68031278,  -11446.00159032, -16047.00169414,
            -20912.46427702, -25898.78333464, -30880.87406771, -35748.39943111, -40403.73628576,
            -44760.97696919, -48745.62840863, -52294.85243412, -55358.20187052, -57898.90076092,
            -59895.83672482, -61346.65187992, -62272.8430536,  -62729.46895957, -62729.46895957,
            -62729.46895957, -62729.46895957, -62729.46895957, -62729.46895957, -62729.46895957,
        },
    .expected_gdab =
        {
            -0.00293509, -0.0034444,  -0.00384158, -0.00414579, -0.00437582,
            -0.0045481,  -0.00467618, -0.00477084, -0.00484045, -0.00489145,
            -0.00492869, -0.0049558,  -0.00497548, -0.00498974, -0.00500006,
            -0.00500751, -0.00501288, -0.00501674, -0.00501952, -0.00502152,
            -0.00502295, -0.00502397, -0.00502471, -0.00502523, -0.00502561,
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
    double r_array_spheres[ARRAY_LEN] = {1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0,
                                         10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0,
                                         19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0};

    double r_array_gdab[ARRAY_LEN] = {
        15.,         18.54166667, 22.08333333, 25.625,      29.16666667, 32.70833333, 36.25,
        39.79166667, 43.33333333, 46.875,      50.41666667, 53.95833333, 57.5,        61.04166667,
        64.58333333, 68.125,      71.66666667, 75.20833333, 78.75,       82.29166667, 85.83333333,
        89.375,      92.91666667, 96.45833333, 100.,
    };

    compute_analytical_spheres(&params_spheres, r_array_spheres, G_analytic_spheres, ARRAY_LEN);
    memcpy(ctx.actual_spheres, G_analytic_spheres, sizeof ctx.actual_spheres);

    compute_analytical_gdab(&params_gdab, r_array_gdab, G_analytic_gdab, ARRAY_LEN);
    memcpy(ctx.actual_gdab, G_analytic_gdab, sizeof ctx.actual_gdab);

    // printf("Analitycal solution:   ");
    // for (size_t i = 0; i < ARRAY_LEN; i++) {
    //     printf("%.8g, ", (G_analytic_gdab[i]));
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
