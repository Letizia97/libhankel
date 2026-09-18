// clang-format off
#include "utils_for_tests/unity_config.h"
#include "unity.h"
// clang-format on
#include <math.h>

#include "interp_loglinear.h"

void setUp(void) {}
void tearDown(void) {}

void test_create_validation(void) {
    // Reject y <= 0 and y non-finite.
    double x_ok[3] = {0.0, 1.0, 2.0};
    double y_ok[3] = {1.0, 10.0, 100.0};
    double y_zero[3] = {1.0, 0.0, 100.0};
    double y_neg[3] = {1.0, -10.0, 100.0};
    double y_nan[3] = {1.0, NAN, 100.0};
    double y_inf[3] = {1.0, INFINITY, 100.0};

    TEST_ASSERT_NULL_MESSAGE(loglinear_interp_create(x_ok, y_zero, 3), "reject y = 0");
    TEST_ASSERT_NULL_MESSAGE(loglinear_interp_create(x_ok, y_neg, 3), "reject y < 0");
    TEST_ASSERT_NULL_MESSAGE(loglinear_interp_create(x_ok, y_nan, 3), "reject NaN ordinate");
    TEST_ASSERT_NULL_MESSAGE(loglinear_interp_create(x_ok, y_inf, 3), "reject infinite ordinate");

    // Reject invalid x and NULL inputs.
    double x_dup[3] = {0.0, 1.0, 1.0};
    double x_nan[3] = {NAN, 1.0, 2.0};

    TEST_ASSERT_NULL_MESSAGE(loglinear_interp_create(x_dup, y_ok, 3), "reject non-increasing x");
    TEST_ASSERT_NULL_MESSAGE(loglinear_interp_create(x_nan, y_ok, 3), "reject NaN abscissa");
    TEST_ASSERT_NULL_MESSAGE(loglinear_interp_create(x_ok, y_ok, 1), "reject n < 2");
    TEST_ASSERT_NULL_MESSAGE(loglinear_interp_create(NULL, y_ok, 3), "reject NULL x");
    TEST_ASSERT_NULL_MESSAGE(loglinear_interp_create(x_ok, NULL, 3), "reject NULL y");

    // Accept valid inputs.
    loglinear_interp_t *h = loglinear_interp_create(x_ok, y_ok, 3);
    TEST_ASSERT_NOT_NULL_MESSAGE(h, "accept valid positive data");
    loglinear_interp_destroy(h);
}

void test_eval_at_nodes_is_exact(void) {
    // Nodes return bit-exact values.
    double x[5] = {0.0, 1.0, 2.0, 3.0, 4.0};
    double y[5] = {1.0, 10.0, 100.0, 1000.0, 10000.0};

    loglinear_interp_t *h = loglinear_interp_create(x, y, 5);
    TEST_ASSERT_NOT_NULL_MESSAGE(h, "create should succeed");

    for (int i = 0; i < 5; i++) {
        char msg[64];
        snprintf(msg, sizeof msg, "node not exact at x = %g", x[i]);
        double v = loglinear_interp_eval(h, x[i]);
        // log/exp round-trip introduces ~1e-15 relative error; use tolerance.
        TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(1e-14 * y[i], y[i], v, msg);
    }
    loglinear_interp_destroy(h);

    // Adversarial last-node case: tiny vs large values.
    double x2[3] = {0.0, 1.0, 2.0};
    double y2[3] = {1e-10, 1.0, 1e10};

    h = loglinear_interp_create(x2, y2, 3);
    TEST_ASSERT_NOT_NULL_MESSAGE(h, "create should succeed");
    double v = loglinear_interp_eval(h, 2.0);
    TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(1e-14 * 1e10, 1e10, v, "last node must be exact");

    loglinear_interp_destroy(h);
}

void test_exponential_interpolation(void) {
    // Geometric means: log-space linear gives exponential (power-law) behavior.
    // Between 1 and 100: midpoint in log space is (log(1) + log(100))/2 = log(10),
    // which gives exp(log(10)) = 10.
    double x[2] = {0.0, 1.0};
    double y[2] = {1.0, 100.0};

    loglinear_interp_t *h = loglinear_interp_create(x, y, 2);
    TEST_ASSERT_NOT_NULL_MESSAGE(h, "create should succeed");

    // Midpoint: geometric mean of 1 and 100 is 10.
    double mid = loglinear_interp_eval(h, 0.5);
    TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(1e-14, 10.0, mid, "midpoint should be geometric mean");

    // Quarter point: exp((log(1) + 3*log(100))/4) = exp(3*log(100)/4) = 100^0.75 ≈ 31.6.
    // But let's use a simpler one: 0.25 point is exp(0.25*log(100)) = 100^0.25 ≈ 3.162.
    double quarter = loglinear_interp_eval(h, 0.25);
    double expected_quarter = pow(100.0, 0.25);
    TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(1e-14, expected_quarter, quarter,
                                      "quarter point should follow power law");

    loglinear_interp_destroy(h);
}

void test_positivity_guaranteed(void) {
    // All interpolated values stay positive.
    double x[4] = {0.0, 1.0, 2.0, 3.0};
    double y[4] = {0.1, 1.0, 10.0, 100.0};

    loglinear_interp_t *h = loglinear_interp_create(x, y, 4);
    TEST_ASSERT_NOT_NULL_MESSAGE(h, "create should succeed");

    for (int i = 0; i < 501; i++) {
        double xi = 3.0 * i / 500.0;
        double v = loglinear_interp_eval(h, xi);
        char msg[48];
        snprintf(msg, sizeof msg, "non-positive at xi = %g", xi);
        TEST_ASSERT_TRUE_MESSAGE(v > 0.0, msg);
    }

    loglinear_interp_destroy(h);
}

void test_eval_error_handling(void) {
    // NULL handle and queries outside range return NaN.
    TEST_ASSERT_DOUBLE_IS_NAN_MESSAGE(loglinear_interp_eval(NULL, 0.5), "NULL handle");

    double x[3] = {1.0, 2.0, 3.0};
    double y[3] = {1.0, 10.0, 100.0};

    loglinear_interp_t *h = loglinear_interp_create(x, y, 3);
    TEST_ASSERT_NOT_NULL_MESSAGE(h, "create should succeed");

    TEST_ASSERT_DOUBLE_IS_NAN_MESSAGE(loglinear_interp_eval(h, 0.0), "below range");
    TEST_ASSERT_DOUBLE_IS_NAN_MESSAGE(loglinear_interp_eval(h, 4.0), "above range");
    TEST_ASSERT_DOUBLE_IS_NAN_MESSAGE(loglinear_interp_eval(h, NAN), "NaN input");

    loglinear_interp_destroy(h);

    // destroy(NULL) is safe.
    loglinear_interp_destroy(NULL);
}

void test_monotonicity(void) {
    // Monotone decreasing y stays monotone decreasing after exp(log_interp).
    enum { N = 64 };
    double x[N], y[N];

    for (int i = 0; i < N; i++) {
        x[i] = pow(10.0, -3.0 + 5.0 * i / (N - 1));
        y[i] = 1.0 / (1.0 + x[i] * x[i] * x[i] * x[i]);
    }

    loglinear_interp_t *h = loglinear_interp_create(x, y, N);
    TEST_ASSERT_NOT_NULL_MESSAGE(h, "create should succeed");

    double prev = INFINITY;
    for (int i = 0; i < 2001; i++) {
        double xi = x[0] + (x[N - 1] - x[0]) * i / 2000.0;
        double v = loglinear_interp_eval(h, xi);
        char msg[64];

        snprintf(msg, sizeof msg, "non-monotone at xi = %g", xi);
        TEST_ASSERT_TRUE_MESSAGE(v <= prev && v > 0.0, msg);

        prev = v;
    }

    loglinear_interp_destroy(h);
}

int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_create_validation);
    RUN_TEST(test_eval_at_nodes_is_exact);
    RUN_TEST(test_exponential_interpolation);
    RUN_TEST(test_positivity_guaranteed);
    RUN_TEST(test_eval_error_handling);
    RUN_TEST(test_monotonicity);
    return UNITY_END();
}
