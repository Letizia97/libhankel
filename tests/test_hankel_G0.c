// clang-format off
#include "utils_for_tests/unity_config.h"
#include "unity.h"
// clang-format on
#include <libhankel.h>
#include <math.h>
#include <stddef.h>

#include "form_factors.h"

form_factor_ctx ctx_spheres;
form_factor_ctx ctx_gdab;

typedef struct {
    double actual_spheres;
    double expected_spheres;
    double actual_gdab;
    double expected_gdab;
} TestContext;

TestContext ctx = {
    .expected_spheres = 394781.721568501,
    .actual_spheres = 0,
    .expected_gdab = 0.0315827340834857,
    .actual_gdab = 0,
};

static double exponential_ff(double q, void *f_ctx) {
    (void)f_ctx;
    return exp(-q);
}

void setUp(void) {
    static double params_spheres[] = {10.0, 1.0};
    ctx_spheres.params = params_spheres;

    static double params_gdab[] = {10.0, 0.5, 1e-4};
    ctx_gdab.params = params_gdab;

    ctx.actual_spheres = hankel_G0(form_factor_sphere, &ctx_spheres, 1e-9);
    ctx.actual_gdab = hankel_G0(form_factor_g_dab, &ctx_gdab, 1e-9);
}

void tearDown(void) {}

void test_hankel_G0_spheres(void) {
    TEST_ASSERT_DOUBLE_WITHIN(1e-3, ctx.expected_spheres, ctx.actual_spheres);
}

void test_hankel_G0_gdab(void) {
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, ctx.expected_gdab, ctx.actual_gdab);
}

void test_hankel_G0_exponential_analytical(void) {
    /* int_0^inf q * exp(-q) dq = Gamma(2) = 1 */
    double result = hankel_G0(exponential_ff, NULL, 1e-9);
    TEST_ASSERT_DOUBLE_WITHIN(1e-9, 1.0, result);
}

int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_hankel_G0_spheres);
    RUN_TEST(test_hankel_G0_gdab);
    RUN_TEST(test_hankel_G0_exponential_analytical);
    return UNITY_END();
}
