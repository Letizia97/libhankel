// clang-format off
#include "utils_for_tests/unity_config.h"
#include "unity.h"
// clang-format on
#include <math.h>
#include <stdio.h>

/*
Tests for the piecewise linear interpolator in src/utils/interp_linear.c.

Deliberately parallel to test_interp_cubic.c, since the two interpolators
publish the same API and the same NaN-outside-range contract. The tests that
have no cubic counterpart are the three guarantees linear makes and PCHIP does
not: exactness at the nodes down to the last bit, containment within the
bracketing pair of ordinates, and monotonicity in xi.
*/
#include "interp_linear.h"

void setUp(void) {}
void tearDown(void) {}

void test_eval_at_node_returns_exact_value(void) {
    /*
    Stronger than the cubic equivalent, which only asks for 1e-12: linear
    interpolation short-circuits the nodes before doing any arithmetic, so the
    tabulated value must come back bit for bit. A tolerance here would let a
    regression in that short-circuit through unnoticed.
    */
    double x[5] = {0.0, 1.0, 2.0, 3.0, 4.0};
    double y[5] = {0.0, 1.0, 4.0, 9.0, 16.0};
    int i;

    linear_interp_t *h = linear_interp_create(x, y, 5);
    TEST_ASSERT_NOT_NULL_MESSAGE(h, "create should succeed for 5 valid increasing points");

    for (i = 0; i < 5; i++) {
        char msg[64];
        snprintf(msg, sizeof msg, "node not reproduced exactly at x = %g", x[i]);
        TEST_ASSERT_TRUE_MESSAGE(linear_interp_eval(h, x[i]) == y[i], msg);
    }

    linear_interp_destroy(h);
}

void test_eval_at_last_node_is_exact(void) {
    /*
    The last node is the one find_interval gets structurally wrong: every
    x[mid] is below it, so lo walks up to n-2 and the returned interval is
    [x[n-2], x[n-1]] with xi on its closed upper end. Interior nodes are safe
    for free -- they bind to the left, giving t = 0 and an exact f0 + 0 -- so
    the top node is the only place the upper guards do any work.

    The ordinates are chosen adversarially, and it took two attempts to get
    them right -- mutation testing found both mistakes. A tidy table like
    {10, 20, 30, 40} lands on the node by arithmetic alone and passes with
    every guard deleted. A wildly split pair like {0.5, 1e-20} misses, but
    misses *past* f1, so blend's clamp pulls it back onto the node and the
    test still cannot fail.

    1.0 -> 1e-16 is the case that discriminates: f1 - f0 rounds to
    -0.9999999999999999, the sum comes back as 2^-53, and that lands strictly
    inside [1e-16, 1.0] where blend's clamp has no reason to touch it. Delete
    the xi == x[i + 1] short-circuit and this assertion reads
    1.11e-16 != 1e-16.
    */
    double x[3] = {0.0, 1.0, 2.0};
    double y[3] = {2.0, 1.0, 1e-16};

    linear_interp_t *h = linear_interp_create(x, y, 3);
    TEST_ASSERT_NOT_NULL_MESSAGE(h, "create should succeed for 3 valid increasing points");

    TEST_ASSERT_TRUE_MESSAGE(linear_interp_eval(h, 2.0) == 1e-16,
                             "the top node must come back exactly, not interpolated");

    linear_interp_destroy(h);
}

void test_interpolates_linearly_between_nodes(void) {
    /*
    The defining behaviour. A midpoint must be the mean of the bracketing
    ordinates, and a quarter point the corresponding weighted mean. Both are
    representable exactly here, so an exact comparison is fair and catches a
    reversed weight -- which a symmetric midpoint-only test would not.
    */
    double x[3] = {0.0, 2.0, 4.0};
    double y[3] = {0.0, 8.0, 16.0};

    linear_interp_t *h = linear_interp_create(x, y, 3);
    TEST_ASSERT_NOT_NULL_MESSAGE(h, "create should succeed for 3 valid increasing points");

    TEST_ASSERT_TRUE_MESSAGE(linear_interp_eval(h, 1.0) == 4.0, "midpoint of [0, 2] should be 4");
    TEST_ASSERT_TRUE_MESSAGE(linear_interp_eval(h, 0.5) == 2.0, "quarter point of [0, 2] should be 2");
    TEST_ASSERT_TRUE_MESSAGE(linear_interp_eval(h, 3.0) == 12.0, "midpoint of [2, 4] should be 12");

    linear_interp_destroy(h);
}

void test_two_point_table_is_accepted(void) {
    /*
    n = 2 is the documented minimum and the degenerate case for the bisection:
    hi - lo is already 1, so the loop body never runs and find_interval must
    return 0 without touching x[mid]. An off-by-one in the loop guard shows up
    here and nowhere else.
    */
    double x[2] = {0.0, 1.0};
    double y[2] = {5.0, 7.0};

    linear_interp_t *h = linear_interp_create(x, y, 2);
    TEST_ASSERT_NOT_NULL_MESSAGE(h, "create should accept the two-point minimum");

    TEST_ASSERT_TRUE_MESSAGE(linear_interp_eval(h, 0.0) == 5.0, "lower node of a two-point table");
    TEST_ASSERT_TRUE_MESSAGE(linear_interp_eval(h, 1.0) == 7.0, "upper node of a two-point table");
    TEST_ASSERT_TRUE_MESSAGE(linear_interp_eval(h, 0.5) == 6.0, "midpoint of a two-point table");

    linear_interp_destroy(h);
}

void test_create_returns_null_for_non_increasing_x(void) {
    /*
    A repeated abscissa would make the interval width zero and the weight a
    division by zero. Unlike the cubic wrapper, which is catching a Boost
    exception, this rejection is the validation loop's own work.
    */
    double x[5] = {0.0, 1.0, 1.0, 3.0, 4.0};
    double y[5] = {0.0, 1.0, 4.0, 9.0, 16.0};

    TEST_ASSERT_NULL_MESSAGE(linear_interp_create(x, y, 5),
                             "create should reject x that is not strictly increasing");
}

void test_create_returns_null_for_too_few_points(void) {
    /*
    One point defines no interval at all. The threshold is 2 rather than the
    cubic's 4, so this also pins the difference between the two headers.
    */
    double x[1] = {0.0};
    double y[1] = {1.0};

    TEST_ASSERT_NULL_MESSAGE(linear_interp_create(x, y, 1), "create should return null if n < 2");
}

void test_create_returns_null_for_null_input(void) {
    /*
    Everything but the NULL is valid, so nothing else can reject the call.
    */
    double x[3] = {0.0, 1.0, 2.0};
    double y[3] = {0.0, 1.0, 4.0};

    TEST_ASSERT_NULL_MESSAGE(linear_interp_create(NULL, y, 3), "create should reject a NULL x pointer");
    TEST_ASSERT_NULL_MESSAGE(linear_interp_create(x, NULL, 3), "create should reject a NULL y pointer");
}

void test_create_returns_null_for_non_finite_x(void) {
    /*
    Two separate guards, and each needs its own case. NaN fails the
    strictly-increasing comparison, so it is caught even at index 0 where
    there is no predecessor to compare against -- but only by the isfinite
    call. An infinity passes the increasing test happily and would make every
    interval width infinite and every weight zero, which is what isfinite is
    there to stop.
    */
    double x_nan[3] = {NAN, 1.0, 2.0};
    double x_inf[3] = {0.0, 1.0, INFINITY};
    double y[3] = {0.0, 1.0, 4.0};

    TEST_ASSERT_NULL_MESSAGE(linear_interp_create(x_nan, y, 3), "create should reject a NaN abscissa");
    TEST_ASSERT_NULL_MESSAGE(linear_interp_create(x_inf, y, 3), "create should reject an infinite abscissa");
}

void test_eval_returns_nan_for_null_handle(void) {
    /*
    Lets a caller who ignored a NULL from create keep going: the bad handle
    poisons the result with NaN instead of dereferencing through it.
    */
    TEST_ASSERT_DOUBLE_IS_NAN_MESSAGE(linear_interp_eval(NULL, 2.0),
                                      "eval should return NaN for NULL handle");
}

void test_eval_returns_nan_outside_range(void) {
    /*
    The contract that makes this unsafe to hand straight to hankel_transform,
    and the reason tabulated_ff exists. Both ends are asserted because they are
    separate comparisons and a broken one could hide behind a working one.
    */
    double x[4] = {1.0, 2.0, 3.0, 4.0};
    double y[4] = {1.0, 4.0, 9.0, 16.0};

    linear_interp_t *h = linear_interp_create(x, y, 4);
    TEST_ASSERT_NOT_NULL_MESSAGE(h, "create should succeed for 4 valid increasing points");

    TEST_ASSERT_DOUBLE_IS_NAN_MESSAGE(linear_interp_eval(h, 0.0),
                                      "eval should return NaN below the bottom of the range");
    TEST_ASSERT_DOUBLE_IS_NAN_MESSAGE(linear_interp_eval(h, 5.0),
                                      "eval should return NaN above the top of the range");

    linear_interp_destroy(h);
}

void test_eval_returns_nan_for_nan_input(void) {
    /*
    The range check is written negated precisely so a NaN xi falls through it.
    Written as (xi < x[0] || xi > x[n-1]) instead, both comparisons would be
    false for NaN and the point would be handed to the bisection, which would
    place it in some arbitrary interval and return a plausible-looking number.
    */
    double x[4] = {1.0, 2.0, 3.0, 4.0};
    double y[4] = {1.0, 4.0, 9.0, 16.0};

    linear_interp_t *h = linear_interp_create(x, y, 4);
    TEST_ASSERT_NOT_NULL_MESSAGE(h, "create should succeed for 4 valid increasing points");

    TEST_ASSERT_DOUBLE_IS_NAN_MESSAGE(linear_interp_eval(h, NAN),
                                      "eval should return NaN for a NaN abscissa");

    linear_interp_destroy(h);
}

void test_destroy_null_is_safe(void) {
    /*
    Pins the contract that lets cleanup paths call destroy without checking h,
    exactly as free() allows. Reaching the next line is the pass.
    */
    linear_interp_destroy(NULL);

    TEST_PASS();
}

void test_no_overshoot_on_step_data(void) {
    /*
    The cubic file's counterpart test allows 1e-12 of slop because PCHIP only
    promises not to ring. Linear promises the result never leaves the
    bracketing pair at all, so this asserts containment exactly. The clamp at
    the end of blend is what makes an exact bound defensible.

    Sampled between the nodes on purpose: every interpolator reproduces its own
    data points, so a test that only looked at the nodes would pass for an
    implementation that overshoots in between.
    */
    double x[6] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    double y[6] = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
    const int nsamples = 501;
    int i;

    linear_interp_t *h = linear_interp_create(x, y, 6);
    TEST_ASSERT_NOT_NULL_MESSAGE(h, "create should succeed for 6 valid increasing points");

    for (i = 0; i < nsamples; i++) {
        /* Recomputed from i rather than accumulated: xi += step would drift
         * off the end of the range and evaluate to NaN on the last sample. */
        double xi = 5.0 * i / (nsamples - 1);
        double v = linear_interp_eval(h, xi);
        char msg[64];

        /* Every iteration reports the same __LINE__, so the sample point has
         * to go in the message or a failure says nothing about where. */
        snprintf(msg, sizeof msg, "value left [0, 1] at xi = %g", xi);

        TEST_ASSERT_TRUE_MESSAGE(v >= 0.0 && v <= 1.0, msg);
    }

    linear_interp_destroy(h);
}

void test_stays_within_bracketing_pair_across_zero(void) {
    /*
    Ordinates that straddle zero are the case blend used to special-case with a
    weighted average, on the grounds that the correction form could overflow.
    That branch is gone -- overflow needs ordinates near 1e308, which no form
    factor reaches -- so this test guards what actually replaced it: the final
    clamp keeping the result inside the bracketing pair even when the two ends
    have opposite signs. A scattering amplitude crossing zero is the realistic
    version of this.
    */
    double x[4] = {0.0, 1.0, 2.0, 3.0};
    double y[4] = {-2.0, -1.0, 1.0, 2.0};
    const int nsamples = 501;
    int i;

    linear_interp_t *h = linear_interp_create(x, y, 4);
    TEST_ASSERT_NOT_NULL_MESSAGE(h, "create should succeed for 4 valid increasing points");

    for (i = 0; i < nsamples; i++) {
        /* Confined to [1, 2], the one interval whose ordinates change sign. */
        double xi = 1.0 + (double)i / (nsamples - 1);
        double v = linear_interp_eval(h, xi);
        char msg[64];

        snprintf(msg, sizeof msg, "value left [-1, 1] at xi = %g", xi);

        TEST_ASSERT_TRUE_MESSAGE(v >= -1.0 && v <= 1.0, msg);
    }

    linear_interp_destroy(h);
}

void test_is_monotone_in_xi(void) {
    /*
    The guarantee that stops a densely sampled curve wobbling, and the reason
    blend uses f0 + t*(f1 - f0) rather than the weighted average: the weighted
    form's two products round independently, so stepping xi up can step the
    result down by an ulp. Monotone data in, monotone samples out, with no
    tolerance -- a single backwards step is a real failure.

    Log-spaced abscissae and a decaying ordinate, so the sampling resembles a
    tabulated form factor rather than a tidy unit grid.
    */
    enum { N = 64 };
    double x[N], y[N];
    const int nsamples = 20001;
    double prev;
    int i;

    for (i = 0; i < N; i++) {
        x[i] = pow(10.0, -3.0 + 5.0 * i / (N - 1));
        y[i] = 1.0 / (1.0 + x[i] * x[i] * x[i] * x[i]);
    }

    linear_interp_t *h = linear_interp_create(x, y, N);
    TEST_ASSERT_NOT_NULL_MESSAGE(h, "create should succeed for a log-spaced table");

    prev = INFINITY;
    for (i = 0; i < nsamples; i++) {
        double xi = x[0] + (x[N - 1] - x[0]) * i / (nsamples - 1);
        double v = linear_interp_eval(h, xi);
        char msg[80];

        snprintf(msg, sizeof msg, "result increased at xi = %g on decreasing data", xi);
        TEST_ASSERT_TRUE_MESSAGE(v <= prev, msg);

        snprintf(msg, sizeof msg, "non-negative data gave a negative result at xi = %g", xi);
        TEST_ASSERT_TRUE_MESSAGE(v >= 0.0, msg);

        prev = v;
    }

    linear_interp_destroy(h);
}

int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_eval_at_node_returns_exact_value);
    RUN_TEST(test_eval_at_last_node_is_exact);
    RUN_TEST(test_interpolates_linearly_between_nodes);
    RUN_TEST(test_two_point_table_is_accepted);
    RUN_TEST(test_create_returns_null_for_non_increasing_x);
    RUN_TEST(test_create_returns_null_for_too_few_points);
    RUN_TEST(test_create_returns_null_for_null_input);
    RUN_TEST(test_create_returns_null_for_non_finite_x);
    RUN_TEST(test_eval_returns_nan_for_null_handle);
    RUN_TEST(test_eval_returns_nan_outside_range);
    RUN_TEST(test_eval_returns_nan_for_nan_input);
    RUN_TEST(test_destroy_null_is_safe);
    RUN_TEST(test_no_overshoot_on_step_data);
    RUN_TEST(test_stays_within_bracketing_pair_across_zero);
    RUN_TEST(test_is_monotone_in_xi);
    return UNITY_END();
}
