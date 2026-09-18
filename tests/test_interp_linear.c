// clang-format off
#include "utils_for_tests/unity_config.h"
#include "unity.h"
// clang-format on
#include <math.h>

#include "interp_linear.h"

void setUp(void) {}
void tearDown(void) {}

void test_eval_at_nodes_is_exact(void) {
    // Linear short-circuits nodes before arithmetic, so bit-exact reproduction.
    double x[5] = {0.0, 1.0, 2.0, 3.0, 4.0};
    double y[5] = {0.0, 1.0, 4.0, 9.0, 16.0};

    linear_interp_t *h = linear_interp_create(x, y, 5);
    TEST_ASSERT_NOT_NULL_MESSAGE(h, "create should succeed for 5 valid points");

    for (int i = 0; i < 5; i++) {
        char msg[64];
        snprintf(msg, sizeof msg, "node not exact at x = %g", x[i]);
        TEST_ASSERT_TRUE_MESSAGE(linear_interp_eval(h, x[i]) == y[i], msg);
    }

    // Adversarial last-node case: the xi == x[i+1] guard is what keeps this exact.
    double x2[3] = {0.0, 1.0, 2.0};
    double y2[3] = {2.0, 1.0, 1e-16};
    linear_interp_destroy(h);

    h = linear_interp_create(x2, y2, 3);
    TEST_ASSERT_NOT_NULL_MESSAGE(h, "create should succeed for 3 valid points");
    TEST_ASSERT_TRUE_MESSAGE(linear_interp_eval(h, 2.0) == 1e-16, "last node must be exact");

    linear_interp_destroy(h);
}

void test_interpolates_linearly_between_nodes(void) {
    // Exact representable midpoints and quarter points catch reversed weights.
    double x[3] = {0.0, 2.0, 4.0};
    double y[3] = {0.0, 8.0, 16.0};

    linear_interp_t *h = linear_interp_create(x, y, 3);
    TEST_ASSERT_NOT_NULL_MESSAGE(h, "create should succeed");

    TEST_ASSERT_TRUE_MESSAGE(linear_interp_eval(h, 1.0) == 4.0, "midpoint of [0, 2]");
    TEST_ASSERT_TRUE_MESSAGE(linear_interp_eval(h, 0.5) == 2.0, "quarter point");
    TEST_ASSERT_TRUE_MESSAGE(linear_interp_eval(h, 3.0) == 12.0, "midpoint of [2, 4]");

    linear_interp_destroy(h);
}

void test_create_validation(void) {
    // n = 2 minimum, bisection off-by-one shows here.
    double x[2] = {0.0, 1.0};
    double y[2] = {5.0, 7.0};

    linear_interp_t *h = linear_interp_create(x, y, 2);
    TEST_ASSERT_NOT_NULL_MESSAGE(h, "accept n=2 minimum");

    TEST_ASSERT_TRUE_MESSAGE(linear_interp_eval(h, 0.5) == 6.0, "two-point midpoint");
    linear_interp_destroy(h);

    // Non-increasing, non-finite, and null inputs are rejected.
    double x_ok[3] = {0.0, 1.0, 2.0};
    double y_ok[3] = {0.0, 1.0, 4.0};
    double x_dup[3] = {0.0, 1.0, 1.0};
    double x_nan[3] = {NAN, 1.0, 2.0};
    double x_inf[3] = {0.0, 1.0, INFINITY};

    TEST_ASSERT_NULL_MESSAGE(linear_interp_create(x_dup, y_ok, 3), "reject non-increasing x");
    TEST_ASSERT_NULL_MESSAGE(linear_interp_create(x_ok, y_ok, 1), "reject n < 2");
    TEST_ASSERT_NULL_MESSAGE(linear_interp_create(NULL, y_ok, 3), "reject NULL x");
    TEST_ASSERT_NULL_MESSAGE(linear_interp_create(x_ok, NULL, 3), "reject NULL y");
    TEST_ASSERT_NULL_MESSAGE(linear_interp_create(x_nan, y_ok, 3), "reject NaN abscissa");
    TEST_ASSERT_NULL_MESSAGE(linear_interp_create(x_inf, y_ok, 3), "reject infinite abscissa");
}

void test_eval_error_handling(void) {
    // NULL handle and queries outside range return NaN.
    TEST_ASSERT_DOUBLE_IS_NAN_MESSAGE(linear_interp_eval(NULL, 2.0), "NULL handle");

    double x[4] = {1.0, 2.0, 3.0, 4.0};
    double y[4] = {1.0, 4.0, 9.0, 16.0};

    linear_interp_t *h = linear_interp_create(x, y, 4);
    TEST_ASSERT_NOT_NULL_MESSAGE(h, "create should succeed");

    TEST_ASSERT_DOUBLE_IS_NAN_MESSAGE(linear_interp_eval(h, 0.0), "below range");
    TEST_ASSERT_DOUBLE_IS_NAN_MESSAGE(linear_interp_eval(h, 5.0), "above range");
    TEST_ASSERT_DOUBLE_IS_NAN_MESSAGE(linear_interp_eval(h, NAN), "NaN input");

    linear_interp_destroy(h);

    // destroy(NULL) is safe.
    linear_interp_destroy(NULL);
}

void test_containment(void) {
    // blend's clamp ensures result never leaves [f0, f1].
    double x[6] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    double y[6] = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};

    linear_interp_t *h = linear_interp_create(x, y, 6);
    TEST_ASSERT_NOT_NULL_MESSAGE(h, "create should succeed");

    for (int i = 0; i < 501; i++) {
        double xi = 5.0 * i / 500.0;
        double v = linear_interp_eval(h, xi);
        char msg[48];
        snprintf(msg, sizeof msg, "left [0, 1] at xi = %g", xi);
        TEST_ASSERT_TRUE_MESSAGE(v >= 0.0 && v <= 1.0, msg);
    }
    linear_interp_destroy(h);

    // Crossing zero: ordinates straddle zero.
    double x2[4] = {0.0, 1.0, 2.0, 3.0};
    double y2[4] = {-2.0, -1.0, 1.0, 2.0};

    h = linear_interp_create(x2, y2, 4);
    TEST_ASSERT_NOT_NULL_MESSAGE(h, "create should succeed");

    for (int i = 0; i < 501; i++) {
        double xi = 1.0 + (double)i / 500.0;
        double v = linear_interp_eval(h, xi);
        char msg[48];
        snprintf(msg, sizeof msg, "left [-1, 1] at xi = %g", xi);
        TEST_ASSERT_TRUE_MESSAGE(v >= -1.0 && v <= 1.0, msg);
    }
    linear_interp_destroy(h);
}

void test_monotonicity(void) {
    // f0 + t*(f1 - f0) ensures monotonicity; weighted average would wobble.
    enum { N = 64 };
    double x[N], y[N];

    for (int i = 0; i < N; i++) {
        x[i] = pow(10.0, -3.0 + 5.0 * i / (N - 1));
        y[i] = 1.0 / (1.0 + x[i] * x[i] * x[i] * x[i]);
    }

    linear_interp_t *h = linear_interp_create(x, y, N);
    TEST_ASSERT_NOT_NULL_MESSAGE(h, "create should succeed");

    double prev = INFINITY;
    for (int i = 0; i < 2001; i++) {
        double xi = x[0] + (x[N - 1] - x[0]) * i / 2000.0;
        double v = linear_interp_eval(h, xi);
        char msg[64];

        snprintf(msg, sizeof msg, "non-monotone at xi = %g", xi);
        TEST_ASSERT_TRUE_MESSAGE(v <= prev && v >= 0.0, msg);

        prev = v;
    }

    linear_interp_destroy(h);
}

int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_eval_at_nodes_is_exact);
    RUN_TEST(test_interpolates_linearly_between_nodes);
    RUN_TEST(test_create_validation);
    RUN_TEST(test_eval_error_handling);
    RUN_TEST(test_containment);
    RUN_TEST(test_monotonicity);
    return UNITY_END();
}
