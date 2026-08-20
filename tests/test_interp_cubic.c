// clang-format off
#include "utils/unity_config.h"
#include "unity.h"
// clang-format on
#include <math.h>
#include <stdio.h>

/*
Tests for the cubic (PCHIP) interpolator wrapper in src/utils/interp_cubic.cpp.

Includes the header. cubic_interp_t is an opaque type that cannot be
usefully redeclared here, and including the header also exercises its
__cplusplus linkage guard from a C translation unit.
*/
#include "../src/utils/interp_cubic.h"

void setUp(void) {}
void tearDown(void) {}

void test_eval_at_node_returns_exact_value(void) {
    /*
    The defining property of an interpolator: at a data point the spline must
    give back that point's y value. Also a canary for an off-by-one in the
    vector copy inside cubic_interp_create.
    */
    double x[5] = {0.0, 1.0, 2.0, 3.0, 4.0};
    double y[5] = {0.0, 1.0, 4.0, 9.0, 16.0};

    cubic_interp_t *h = cubic_interp_create(x, y, 5);
    TEST_ASSERT_NOT_NULL_MESSAGE(h, "create should succeed for 5 valid increasing points");

    TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(1e-12, 4.0, cubic_interp_eval(h, 2.0),
                                      "spline must pass through the data point at x = 2");

    cubic_interp_destroy(h);
}

void test_create_returns_null_for_non_increasing_x(void) {
    /*
    Boost throws on a repeated abscissa -- two points at the same x give the
    spline an infinite slope. The wrapper's job is to turn that exception into
    a NULL, since a C caller cannot catch it.
    */
    double x[5] = {0.0, 1.0, 1.0, 3.0, 4.0};
    double y[5] = {0.0, 1.0, 4.0, 9.0, 16.0};

    cubic_interp_t *h = cubic_interp_create(x, y, 5);
    TEST_ASSERT_NULL_MESSAGE(h, "create should reject x that is not strictly increasing");

    cubic_interp_destroy(h);
}

void test_create_returns_null_for_too_few_points(void) {
    /*
    A cubic through fewer than four points is underdetermined, so Boost
    refuses to build one. Same conversion as above: exception becomes NULL.
    x is strictly increasing so the point count is the only thing wrong.
    */
    double x[3] = {0.0, 1.0, 2.0};
    double y[3] = {0.0, 1.0, 4.0};

    cubic_interp_t *h = cubic_interp_create(x, y, 3);
    TEST_ASSERT_NULL_MESSAGE(h, "create should return null if n<4");

    cubic_interp_destroy(h);
}

void test_create_returns_null_for_null_input(void) {
    /*
    The one check the wrapper must make itself: a null deref is a segfault,
    not an exception, so the try/catch could never turn it into a NULL.
    Everything but the NULL is valid, so nothing else can reject the call.
    */
    double x[5] = {0.0, 1.0, 2.0, 3.0, 4.0};
    double y[5] = {0.0, 1.0, 4.0, 9.0, 16.0};

    TEST_ASSERT_NULL_MESSAGE(cubic_interp_create(NULL, y, 5),
                             "create should reject a NULL x pointer");
    TEST_ASSERT_NULL_MESSAGE(cubic_interp_create(x, NULL, 5),
                             "create should reject a NULL y pointer");
}

void test_eval_returns_nan_for_null_handle(void) {
    /*
    Guards the same segfault on the eval side, and lets a caller who ignored
    a NULL from create keep going: the bad handle poisons the result with NaN
    instead of killing the process. xi is arbitrary -- it is never read.
    */
    double xi = 2.0;

    TEST_ASSERT_DOUBLE_IS_NAN_MESSAGE(cubic_interp_eval(NULL, xi),
                                      "eval should return NaN for NULL handle");
}

void test_destroy_null_is_safe(void) {
    /*
    Pins the contract that lets cleanup paths call destroy without first
    checking h, exactly as free() allows. No assertion is possible -- destroy
    returns void -- so reaching the next line is the pass, and a crash would
    take the process down before Unity could record a result.
    */
    cubic_interp_destroy(NULL);

    TEST_PASS();
}

void test_eval_returns_nan_outside_range(void) {
    /*
    Boost refuses to extrapolate and throws, so this NaN comes from the catch
    in cubic_interp_eval rather than from its null guard -- which is why the
    handle is asserted non-NULL first. Both ends are separate comparisons
    inside Boost, so a broken one could hide behind a working one.
    */
    double x[5] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double y[5] = {1.0, 4.0, 9.0, 16.0, 25.0};

    cubic_interp_t *h = cubic_interp_create(x, y, 5);
    TEST_ASSERT_NOT_NULL_MESSAGE(h, "create should succeed for 5 valid increasing points");

    TEST_ASSERT_DOUBLE_IS_NAN_MESSAGE(cubic_interp_eval(h, 0.0),
                                      "eval should return NaN below the bottom of the range");
    TEST_ASSERT_DOUBLE_IS_NAN_MESSAGE(cubic_interp_eval(h, 6.0),
                                      "eval should return NaN above the top of the range");

    cubic_interp_destroy(h);
}

void test_eval_at_range_endpoints_is_finite(void) {
    /*
    The mirror of the test above: the range check is inclusive, so x[0] and
    x[n-1] are the one pair of points sitting on the boundary between the NaN
    behaviour and the normal behaviour. Turning a < into a <= inside Boost
    breaks exactly one of these two tests, which is why both exist.
    */
    double x[5] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double y[5] = {1.0, 4.0, 9.0, 16.0, 25.0};

    cubic_interp_t *h = cubic_interp_create(x, y, 5);
    TEST_ASSERT_NOT_NULL_MESSAGE(h, "create should succeed for 5 valid increasing points");

    TEST_ASSERT_DOUBLE_IS_DETERMINATE_MESSAGE(cubic_interp_eval(h, 1.0),
                                              "eval must be finite at the bottom of the range");
    TEST_ASSERT_DOUBLE_IS_DETERMINATE_MESSAGE(cubic_interp_eval(h, 5.0),
                                              "eval must be finite at the top of the range");

    cubic_interp_destroy(h);
}

void test_no_overshoot_on_step_data(void) {
    /*
    The test that justifies choosing PCHIP. Fitted to a step, a C2 spline 
    rings and dips below zero, which for a scattering curve means a
    negative intensity. PCHIP is monotonicity-preserving, so it 
    cannot leave [0, 1] here.

    Sampled between the nodes on purpose: every interpolator reproduces its
    own data points, so a test that only looked at x = 0..5 would pass for a
    spline that overshoots badly in between.
    */
    double x[6] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    double y[6] = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};

    /* Loose enough to absorb rounding in the evaluation, far tighter than the
     * few percent a ringing spline would overshoot by. */
    const double eps = 1e-12;
    const int nsamples = 101;
    int i;

    cubic_interp_t *h = cubic_interp_create(x, y, 6);
    TEST_ASSERT_NOT_NULL_MESSAGE(h, "create should succeed for 6 valid increasing points");

    for (i = 0; i < nsamples; i++) {
        /* Recomputed from i rather than accumulated: xi += step would drift
         * off the end of the range and evaluate to NaN on the last sample. */
        double xi = 5.0 * i / (nsamples - 1);
        double v = cubic_interp_eval(h, xi);
        char msg[64];

        /* Every iteration reports the same __LINE__, so the sample point has
         * to go in the message or a failure says nothing about where. */
        snprintf(msg, sizeof msg, "value left [0, 1] at xi = %g", xi);

        TEST_ASSERT_GREATER_OR_EQUAL_DOUBLE_MESSAGE(0.0 - eps, v, msg);
        TEST_ASSERT_LESS_OR_EQUAL_DOUBLE_MESSAGE(1.0 + eps, v, msg);
    }

    cubic_interp_destroy(h);
}

int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_eval_at_node_returns_exact_value);
    RUN_TEST(test_create_returns_null_for_non_increasing_x);
    RUN_TEST(test_create_returns_null_for_too_few_points);
    RUN_TEST(test_create_returns_null_for_null_input);
    RUN_TEST(test_eval_returns_nan_for_null_handle);
    RUN_TEST(test_destroy_null_is_safe);
    RUN_TEST(test_eval_returns_nan_outside_range);
    RUN_TEST(test_eval_at_range_endpoints_is_finite);
    RUN_TEST(test_no_overshoot_on_step_data);
    return UNITY_END();
}
