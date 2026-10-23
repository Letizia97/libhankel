// clang-format off
#include "utils_for_tests/unity_config.h"
#include "unity.h"
// clang-format on
#include <math.h>
#include <stddef.h>
#include <stdio.h>

#include "form_factors.h"
#include "libhankel.h"
#include "utils_for_tests/test_utils.h"

#define ARRAY_LEN 10

int nu;

form_factor_ctx ctx_spheres;
form_factor_ctx ctx_gdab;

typedef struct {
    double actual_spheres[ARRAY_LEN];
    double expected_spheres[ARRAY_LEN];
    double actual_gdab[ARRAY_LEN];
    double expected_gdab[ARRAY_LEN];
} TestContext;

TestContext ctx = {
    .expected_spheres = {
        387507.04059582,
        371194.98373392,
        349050.87415468,
        322866.86846381,
        293957.90994016,
        263387.30909878,
        232057.35852631,
        200753.92429196,
        170170.33481168,
        140920.02447442,
    },
    .actual_spheres = {0},
    .expected_gdab = {
        0.013140997753292,
        0.0099409586092683,
        0.0074453617175494,
        0.0055339563557995,
        0.0040886406567342,
        0.0030061670136962,
        0.0022014299007379,
        0.0016066904958907,
        0.0011692607977510,
        0.00084881608987727,
    },
    .actual_gdab = {0},
};

void setUp(void) {
    nu = 0;

    static double params_spheres[] = {10.0, 1.0};
    ctx_spheres.params = params_spheres;

    static double params_gdab[] = {10.0, 0.5, 1e-4};
    ctx_gdab.params = params_gdab;

    double x_spheres[ARRAY_LEN] = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};

    double x_gdab[ARRAY_LEN] = {
        15., 18.54166667, 22.08333333, 25.625,      29.16666667,
        32.70833333, 36.25, 39.79166667, 43.33333333, 46.875,
    };

    hankel_transform_fftlog(nu, form_factor_sphere, x_spheres, ARRAY_LEN,
                            (void *)&ctx_spheres, ctx.actual_spheres);

    hankel_transform_fftlog(nu, form_factor_g_dab, x_gdab, ARRAY_LEN,
                            (void *)&ctx_gdab, ctx.actual_gdab);
}

void tearDown(void) {}

void test_hankel_FFTLog_regression_spheres(void) {
    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        TEST_ASSERT_DOUBLE_WITHIN(1e-4, ctx.expected_spheres[i], ctx.actual_spheres[i]);
    }
}

void test_hankel_FFTLog_regression_gdab(void) {
    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        TEST_ASSERT_DOUBLE_WITHIN(1e-8, ctx.expected_gdab[i], ctx.actual_gdab[i]);
    }
}

int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_hankel_FFTLog_regression_spheres);
    RUN_TEST(test_hankel_FFTLog_regression_gdab);
    return UNITY_END();
}
