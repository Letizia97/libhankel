// clang-format off
#include "utils_for_tests/unity_config.h"
#include "unity.h"
// clang-format on
#include <math.h>
#include <stdio.h>

/*
Tests for the tabulated form factor in src/utils/tabulated_ff.c.

The type exists to make a spline safe to hand to hankel_transform, so most of
what is worth testing is refusal: every table rejected here would otherwise
produce NaN or a plausible wrong answer under a status code of 0.
*/
#include "tabulated_ff.h"

/* Long enough that the fit window (a quarter of the table, capped at 40) is
 * the cap rather than the quarter, which is the case real tables hit. */
#define N_TABLE 200

/* Exponent of the default table: Porod's value for a sharp interface, and
 * comfortably above the 3/2 convergence threshold. */
#define TABLE_P 4.0

#define Q_FIRST 1.0e-3
#define Q_LAST 1.0e1

static double table_q[N_TABLE];
static double table_f[N_TABLE];

/* A pure power law, log-spaced. Pure so the fitted exponent has an exact known
 * answer: anything but TABLE_P is the fit being wrong, not the data curving. */
void setUp(void) {
    for (int i = 0; i < N_TABLE; i++) {
        table_q[i] = Q_FIRST * pow(Q_LAST / Q_FIRST, (double)i / (N_TABLE - 1));
        table_f[i] = pow(table_q[i], -TABLE_P);
    }
}

void tearDown(void) {}

/* Handle from the default table, failing the test if create does not succeed.
 * Every evaluation test below starts this way. */
static tabulated_ff_t *make(tabulated_tail_t tail, double exponent) {
    tabulated_ff_t *t = NULL;
    TEST_ASSERT_EQUAL_INT(0, tabulated_ff_create(table_q, table_f, N_TABLE, tail, exponent, &t));
    return t;
}

/* Status from a create that is expected to fail: out stays NULL on every
 * failure path, so there is nothing to destroy. Rejection tests only. */
static int reject(const double *q, const double *f) {
    tabulated_ff_t *t = NULL;
    return tabulated_ff_create(q, f, 5, TABULATED_TAIL_ZERO, 0.0, &t);
}

/* The fitted exponent is not exposed, so it is read back off the tail: at twice
 * q_hi a power law has fallen by 2^-p, one sixteenth for p = 4. That ratio is
 * the only observable the fit produces, and exactly what the transform depends
 * on. q_hi comes off the table, not from Q_LAST: log-spacing leaves them a
 * couple of ULPs apart, enough to land on the wrong side of the join. */
static double tail_drop_at_double(tabulated_ff_t *t) {
    const double q_hi = table_q[N_TABLE - 1];
    return tabulated_ff_eval(2.0 * q_hi, t) / tabulated_ff_eval(q_hi, t);
}

/* --------------------------------------------------------------------------
 * Rejecting tables that cannot be used
 * ----------------------------------------------------------------------- */

void test_create_rejects_null_arguments(void) {
    /* A null deref is a segfault rather than a bad answer, so these can only be
     * caught here. out is checked too, since it is written before anything
     * else. */
    tabulated_ff_t *t = NULL;

    TEST_ASSERT_EQUAL_INT_MESSAGE(
        -13, tabulated_ff_create(NULL, table_f, N_TABLE, TABULATED_TAIL_ZERO, 0.0, &t),
        "NULL q must be rejected");
    TEST_ASSERT_EQUAL_INT_MESSAGE(
        -13, tabulated_ff_create(table_q, NULL, N_TABLE, TABULATED_TAIL_ZERO, 0.0, &t),
        "NULL f must be rejected");
    TEST_ASSERT_EQUAL_INT_MESSAGE(
        -13, tabulated_ff_create(table_q, table_f, N_TABLE, TABULATED_TAIL_ZERO, 0.0, NULL),
        "NULL out must be rejected");
}

void test_create_leaves_out_null_on_failure(void) {
    /* The contract that lets a caller destroy unconditionally after a failed
     * create. A failure path that returned without clearing out would leave the
     * caller freeing whatever was in that variable before the call. */
    tabulated_ff_t *t = (tabulated_ff_t *)0xDEADBEEF;

    TEST_ASSERT_EQUAL_INT(-13,
                          tabulated_ff_create(table_q, table_f, 2, TABULATED_TAIL_ZERO, 0.0, &t));
    TEST_ASSERT_NULL_MESSAGE(t, "out must be cleared even when create fails");
}

void test_create_rejects_too_few_points(void) {
    /* A cubic through fewer than four points is underdetermined. Caught here
     * rather than in Boost, so the caller gets -13 and a point-count message
     * instead of the -3 that a NULL from the spline would become. */
    tabulated_ff_t *t = NULL;

    TEST_ASSERT_EQUAL_INT(-13,
                          tabulated_ff_create(table_q, table_f, 3, TABULATED_TAIL_ZERO, 0.0, &t));
}

void test_create_rejects_non_increasing_q(void) {
    /* Two points at the same abscissa give a segment of zero width and infinite
     * slope. The ordering is also what the spline's search and eval's range
     * comparisons rely on. */
    double q[5] = {1.0, 2.0, 2.0, 4.0, 5.0};
    double f[5] = {1.0, 0.5, 0.4, 0.2, 0.1};

    TEST_ASSERT_EQUAL_INT(-13, reject(q, f));
}

void test_create_rejects_non_positive_q(void) {
    /* log(q) is taken when fitting the tail, and q <= 0 has no meaning anyway.
     * The implementation only checks q[0], relying on strict increase for the
     * rest, so this is what catches that shortcut breaking. */
    double q[5] = {0.0, 1.0, 2.0, 3.0, 4.0};
    double f[5] = {1.0, 0.5, 0.4, 0.2, 0.1};

    TEST_ASSERT_EQUAL_INT(-13, reject(q, f));
}

void test_create_rejects_nan_in_q(void) {
    /* Caught by the strictly-increasing test rather than an explicit check,
     * since any comparison against NaN is false. Its own test because that is a
     * property of how the comparison is spelled: the more natural
     * `q[j] <= q[j - 1]` silently stops catching this. */
    double q[5] = {1.0, 2.0, NAN, 4.0, 5.0};
    double f[5] = {1.0, 0.5, 0.4, 0.2, 0.1};

    TEST_ASSERT_EQUAL_INT(-13, reject(q, f));
}

void test_create_rejects_infinite_q(void) {
    /* The NaN trick does not extend to infinity - INFINITY > 4.0 is simply true
     * - so this needs a check of its own. Without it q_hi is infinite, the tail
     * is unreachable, and every evaluation is NaN under a status code of 0. */
    double q[5] = {1.0, 2.0, 3.0, 4.0, INFINITY};
    double f[5] = {1.0, 0.5, 0.4, 0.2, 0.1};

    TEST_ASSERT_EQUAL_INT(-13, reject(q, f));
}

void test_create_rejects_non_finite_f(void) {
    /* f is only ever read and never compared, so it gets no free NaN check the
     * way q does. Boost builds a spline through either of these quite happily
     * and then returns NaN from every evaluation that touches the segment. */
    double q[5] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double f_nan[5] = {1.0, 0.5, NAN, 0.2, 0.1};
    double f_inf[5] = {1.0, 0.5, 0.4, INFINITY, 0.1};

    TEST_ASSERT_EQUAL_INT_MESSAGE(-13, reject(q, f_nan), "NaN f rejected");
    TEST_ASSERT_EQUAL_INT_MESSAGE(-13, reject(q, f_inf), "infinite f rejected");
}

void test_create_rejects_unknown_tail(void) {
    /* The enum arrives as an integer from the Python binding, so an
     * out-of-range value is reachable from outside C. Unchecked it would fall
     * through eval's `== TABULATED_TAIL_ZERO` test and behave as a power law. */
    tabulated_ff_t *t = NULL;

    TEST_ASSERT_EQUAL_INT(
        -13, tabulated_ff_create(table_q, table_f, N_TABLE, (tabulated_tail_t)99, 4.0, &t));
}

/* --------------------------------------------------------------------------
 * Rejecting tails that would not converge
 * ----------------------------------------------------------------------- */

void test_convergence_threshold_is_exclusive(void) {
    /* The integrand goes like q^(1/2 - p), and at p = 3/2 exactly it goes like
     * q^-1 and diverges logarithmically - so the threshold has to be exclusive.
     * The pair pins which side each value falls on; a > that became a >= would
     * break one of them. NaN rides along because it is the other way a supplied
     * exponent can be unusable, and it is rejected by the same comparison.
     *
     * The reject cases come first: create clears *out, so a failure after the
     * success below would drop the handle it had just written. */
    tabulated_ff_t *t = NULL;

    TEST_ASSERT_EQUAL_INT_MESSAGE(
        -14, tabulated_ff_create(table_q, table_f, N_TABLE, TABULATED_TAIL_POWER_LAW, 1.5, &t),
        "p = 3/2 itself must be rejected");
    TEST_ASSERT_EQUAL_INT_MESSAGE(
        -14, tabulated_ff_create(table_q, table_f, N_TABLE, TABULATED_TAIL_POWER_LAW, NAN, &t),
        "NaN p must not slip through the comparison");
    TEST_ASSERT_EQUAL_INT_MESSAGE(
        0, tabulated_ff_create(table_q, table_f, N_TABLE, TABULATED_TAIL_POWER_LAW, 1.51, &t),
        "just above 3/2 must be accepted");
    tabulated_ff_destroy(t);
}

void test_zero_tail_ignores_the_exponent(void) {
    /* A truncated tail has no exponent, so the argument is documented as unused.
     * Passing one that would be rejected outright as a power law is the
     * sharpest way to show it is not consulted. */
    tabulated_ff_t *t = make(TABULATED_TAIL_ZERO, 1.0);

    TEST_ASSERT_NOT_NULL(t);
    tabulated_ff_destroy(t);
}

void test_create_rejects_exponent_fitted_from_a_truncated_table(void) {
    /* The most valuable test here. A table that stops before the curve reaches
     * its asymptote has a shallow local slope, and the fit can only report the
     * slope it was shown.
     *
     * A Debye-Anderson-Brumberger curve, true asymptote p = 4, sampled only out
     * to q * xi = 0.5, so the fit sees well under 1. Without the guard every
     * transform using it would return plausible, meaningless numbers. */
    const double xi = 10.0;
    double q[N_TABLE];
    double f[N_TABLE];
    tabulated_ff_t *t = NULL;

    for (int i = 0; i < N_TABLE; i++) {
        q[i] = 1.0e-3 * pow(0.05 / 1.0e-3, (double)i / (N_TABLE - 1));
        double u = q[i] * xi;
        f[i] = 1.0 / ((1.0 + u * u) * (1.0 + u * u));
    }

    TEST_ASSERT_EQUAL_INT_MESSAGE(
        -14, tabulated_ff_create(q, f, N_TABLE, TABULATED_TAIL_POWER_LAW, 0.0, &t),
        "a table truncated before the asymptote must not fit a usable exponent");
}

/* --------------------------------------------------------------------------
 * The fit itself
 * ----------------------------------------------------------------------- */

void test_fitted_exponent_recovers_a_known_power_law(void) {
    /* Tight, because the table is an exact power law: the only error is
     * rounding in the sums, not curvature the fit had to average over. */
    tabulated_ff_t *t = make(TABULATED_TAIL_POWER_LAW, 0.0);

    TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(1e-9, pow(2.0, -TABLE_P), tail_drop_at_double(t),
                                      "fitted exponent must reproduce the table's own power law");

    tabulated_ff_destroy(t);
}

void test_fit_skips_non_positive_values(void) {
    /* Background-subtracted data goes negative at high q, precisely where the
     * fit window sits, and log() of a negative is NaN - which would propagate
     * through all four running sums. Only the last three points are spoiled,
     * leaving well above the three survivors needed, so the exponent must come
     * back unchanged. */
    table_f[N_TABLE - 1] = -1.0e-5;
    table_f[N_TABLE - 2] = -2.0e-5;
    table_f[N_TABLE - 3] = 0.0;

    tabulated_ff_t *t = make(TABULATED_TAIL_POWER_LAW, 0.0);

    TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(1e-9, pow(2.0, -TABLE_P), tail_drop_at_double(t),
                                      "the surviving points still describe the same power law");

    tabulated_ff_destroy(t);
}

/* --------------------------------------------------------------------------
 * Evaluation
 * ----------------------------------------------------------------------- */

void test_eval_reproduces_the_tabulated_values(void) {
    /* The defining property of an interpolator, and a canary for an off-by-one
     * in the copy into the spline. Both endpoints are included because eval
     * answers those from its cached values, a different branch from the
     * interior point. */
    tabulated_ff_t *t = make(TABULATED_TAIL_POWER_LAW, 0.0);

    TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(1e-12 * table_f[0], table_f[0],
                                      tabulated_ff_eval(table_q[0], t), "at the first point");
    TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(1e-12 * table_f[100], table_f[100],
                                      tabulated_ff_eval(table_q[100], t), "at an interior point");
    TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(1e-12 * table_f[N_TABLE - 1], table_f[N_TABLE - 1],
                                      tabulated_ff_eval(table_q[N_TABLE - 1], t),
                                      "at the last point");

    tabulated_ff_destroy(t);
}

void test_eval_below_the_table_is_flat(void) {
    /* The low-q plateau, which is not configurable. Two very different q values
     * both far below the table must agree, since a rule that accidentally
     * extrapolated would separate them. */
    tabulated_ff_t *t = make(TABULATED_TAIL_POWER_LAW, 0.0);

    TEST_ASSERT_EQUAL_DOUBLE(table_f[0], tabulated_ff_eval(Q_FIRST * 0.5, t));
    TEST_ASSERT_EQUAL_DOUBLE(table_f[0], tabulated_ff_eval(1.0e-30, t));

    tabulated_ff_destroy(t);
}

void test_eval_above_the_table_follows_the_power_law(void) {
    /* Three decades out, where an error in the exponent is unmissable. Supplied
     * rather than fitted, so the expected value is arithmetic and not a second
     * measurement of the fit. */
    tabulated_ff_t *t = make(TABULATED_TAIL_POWER_LAW, 3.0);

    double expected = table_f[N_TABLE - 1] * pow(1000.0, -3.0);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12 * expected, expected,
                              tabulated_ff_eval(1000.0 * table_q[N_TABLE - 1], t));

    tabulated_ff_destroy(t);
}

void test_eval_above_the_table_is_zero_for_a_truncated_tail(void) {
    /* Exactly zero, not merely small: a tail that leaked a residual value would
     * still put the wrong thing into the integrand. */
    tabulated_ff_t *t = make(TABULATED_TAIL_ZERO, 0.0);
    const double q_hi = table_q[N_TABLE - 1];

    TEST_ASSERT_EQUAL_DOUBLE(0.0, tabulated_ff_eval(q_hi * 1.001, t));
    TEST_ASSERT_EQUAL_DOUBLE(0.0, tabulated_ff_eval(1.0e20, t));

    /* Still the tabulated value at the join, since the comparison there is
     * inclusive - the cliff starts above q_hi, not at it. With a cliff, the few
     * ULPs between q_hi and Q_LAST are the difference between f_hi and zero. */
    TEST_ASSERT_EQUAL_DOUBLE(table_f[N_TABLE - 1], tabulated_ff_eval(q_hi, t));

    tabulated_ff_destroy(t);
}

void test_power_law_tail_joins_continuously(void) {
    /* What the ratio form of the tail buys: at q_hi the ratio is exactly 1, so
     * the tail starts on the spline's last value with no step. A * pow(q, -p)
     * would leave a jump wherever A was rounded, and a step in the integrand
     * becomes ringing across the whole transform. */
    tabulated_ff_t *t = make(TABULATED_TAIL_POWER_LAW, 0.0);

    /* Read off the table, not from Q_LAST: a few ULPs would put both samples in
     * the tail and pass without ever touching the join this claims to test. */
    const double q_hi = table_q[N_TABLE - 1];
    double at_join = tabulated_ff_eval(q_hi, t);
    double just_above = tabulated_ff_eval(q_hi * (1.0 + 1.0e-10), t);

    TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(1e-8 * at_join, at_join, just_above,
                                      "the tail must not step away from the last tabulated value");

    tabulated_ff_destroy(t);
}

void test_eval_is_finite_across_the_whole_strategy_range(void) {
    /* The contract that makes this type usable at all. DHT_Anderson_801 asks
     * for q from about 1e-15 to 3e+20, and a single NaN anywhere in that span
     * poisons the transform while hankel_transform still reports success. */
    tabulated_ff_t *t = make(TABULATED_TAIL_POWER_LAW, 0.0);

    for (int i = -20; i <= 20; i++) {
        double q = pow(10.0, (double)i);
        char msg[64];

        /* Every iteration reports the same __LINE__, so the exponent has to go
         * in the message or a failure says nothing about where. */
        snprintf(msg, sizeof msg, "not finite at q = 1e%d", i);
        TEST_ASSERT_DOUBLE_IS_DETERMINATE_MESSAGE(tabulated_ff_eval(q, t), msg);
    }

    tabulated_ff_destroy(t);
}

void test_eval_returns_nan_for_null_handle(void) {
    /* Lets a caller who ignored a failed create keep going: the bad handle
     * poisons the result with NaN instead of killing the process. */
    TEST_ASSERT_DOUBLE_IS_NAN(tabulated_ff_eval(1.0, NULL));
}

void test_destroy_null_is_safe(void) {
    /* Pins the contract that lets cleanup paths call destroy without first
     * checking, exactly as free() allows. No assertion is possible - destroy
     * returns void - so reaching the next line is the pass. */
    tabulated_ff_destroy(NULL);

    TEST_PASS();
}

int main(void) {
    UNITY_BEGIN();

    RUN_TEST(test_create_rejects_null_arguments);
    RUN_TEST(test_create_leaves_out_null_on_failure);
    RUN_TEST(test_create_rejects_too_few_points);
    RUN_TEST(test_create_rejects_non_increasing_q);
    RUN_TEST(test_create_rejects_non_positive_q);
    RUN_TEST(test_create_rejects_nan_in_q);
    RUN_TEST(test_create_rejects_infinite_q);
    RUN_TEST(test_create_rejects_non_finite_f);
    RUN_TEST(test_create_rejects_unknown_tail);

    RUN_TEST(test_convergence_threshold_is_exclusive);
    RUN_TEST(test_zero_tail_ignores_the_exponent);
    RUN_TEST(test_create_rejects_exponent_fitted_from_a_truncated_table);

    RUN_TEST(test_fitted_exponent_recovers_a_known_power_law);
    RUN_TEST(test_fit_skips_non_positive_values);

    RUN_TEST(test_eval_reproduces_the_tabulated_values);
    RUN_TEST(test_eval_below_the_table_is_flat);
    RUN_TEST(test_eval_above_the_table_follows_the_power_law);
    RUN_TEST(test_eval_above_the_table_is_zero_for_a_truncated_tail);
    RUN_TEST(test_power_law_tail_joins_continuously);
    RUN_TEST(test_eval_is_finite_across_the_whole_strategy_range);
    RUN_TEST(test_eval_returns_nan_for_null_handle);

    RUN_TEST(test_destroy_null_is_safe);

    return UNITY_END();
}
