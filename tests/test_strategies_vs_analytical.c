/*
Accuracy of every Hankel strategy measured against the analytical GDAB
correlation function.

The other test files pin each strategy to numbers a previous run produced, so
they catch a change but cannot say whether the numbers were right to begin
with. Here the reference is @ref compute_analytical_gdab, a closed form, so a
failure means the strategy really is inaccurate.

Two conventions have to be reconciled before the two sides are comparable:

  * hankel_transform returns T(z) = int f(q) J0(qz) q dq. The correlation
    function carries a 1/2pi that the library deliberately leaves to the
    caller, so G(z) = T(z) / 2pi.
  * compute_analytical_gdab returns eta^2 [G(z) - G(0)], the quantity a SESANS
    measurement is sensitive to. G(0) is a z-independent constant that no
    quadrature can know about, so the test subtracts it from the numerical
    side.

Both parameter sets are checked over a z range wide enough to reach the tail,
and one of them uses a Hurst exponent for which nu = H + 1/2 is not an
integer - the case that has to go through a Bessel function of real order.
*/

// clang-format off
#include "utils_for_tests/unity_config.h"
#include "unity.h"
// clang-format on
#include <math.h>
#include <stddef.h>
#include <stdio.h>

#include "form_factors.h"
#include "libhankel.h"
#include "src/utils/analytical_form_factors.h"
#include "utils_for_tests/test_utils.h"

#define MAX_LEN 19

/* One strategy and the accuracy it is expected to reach on this integrand.
 *
 * The tolerances are relative, and were set by measuring the largest error
 * each strategy produced across all three cases below and rounding up to the
 * next power of ten. That leaves uneven headroom - about 5x for
 * Fixed_DE_Ogata and 8x for DHT_Key_51, but 300x for QWE_Chave, whose measured
 * error is already near machine precision. The tight ones are tight on
 * purpose: rounding them up another decade would stop them detecting anything.
 *
 * Two things this does not do. It does not assert that a strategy is as
 * accurate as it ought to be, only that it is no worse than it was when the
 * numbers were taken (gcc 11.4, Boost from the system), so a different libm
 * could in principle shift the tightest entries. And it is not a ranking of
 * the strategies in general: the GDAB form factor is smooth and monotonically
 * decaying, the easiest case for a digital filter. */
typedef struct {
    const char *name;
    double rel_tol;
} strategy_case;

static const strategy_case STRATEGIES[] = {
    {"DHT_Guptasarma", 1e-10},   {"DHT_Guptasarma_Fast", 1e-8},
    {"DHT_Key_51", 1e-6},        {"DHT_Key_101", 1e-9},
    {"DHT_Key_201", 1e-13},      {"DHT_Anderson_801", 1e-7},
    {"Fixed_DE_Ogata", 1e-2},    {"Adaptive_DE_Ooura", 1e-13},
    {"QWE_Key", 1e-10},          {"QWE_Chave", 1e-10},
};

#define N_STRATEGIES (sizeof(STRATEGIES) / sizeof(STRATEGIES[0]))

/* A set of GDAB parameters together with the abscissae to test it on. */
typedef struct {
    const char *label;
    double A;
    double H;
    double eta;
    size_t len;
    double z[MAX_LEN];
} gdab_case;

/* nu = H + 1/2 = 1, an integer. Same parameters and grid as the regression
 * tests, so the two files can be compared point by point. */
static const gdab_case CASE_INTEGER_ORDER = {
    .label = "H=0.5 (integer nu)",
    .A = 10.0,
    .H = 0.5,
    .eta = 1e-4,
    .len = 19,
    .z = {15., 18.54166667, 22.08333333, 25.625, 29.16666667, 32.70833333, 36.25, 39.79166667,
          43.33333333, 46.875, 50.41666667, 53.95833333, 57.5, 61.04166667, 64.58333333, 68.125,
          71.66666667, 75.20833333, 78.75},
};

/* nu = H + 1/2 = 0.8, not an integer, and z reaching u = z/A = 50. Boost's
 * cyl_bessel_k covers real orders, but the Kummer-function detour that used to
 * stand in for it here lost every significant digit past u ~ 15, so the tail
 * of this case is the point of it. */
static const gdab_case CASE_REAL_ORDER = {
    .label = "H=0.3 (real nu), tail to u=50",
    .A = 10.0,
    .H = 0.3,
    .eta = 1e-4,
    .len = 12,
    .z = {10., 25., 50., 80., 120., 160., 200., 260., 320., 380., 440., 500.},
};

/* H just above its lower bound of -1/2, where the integrand q*I(q) decays as
 * q^-(2+2H) = q^-1.02 and is therefore only barely convergent. This is the
 * hardest tail any strategy has to cover, and the relative error of the
 * transform itself degrades steadily towards it: Adaptive_DE_Ooura loses about
 * three orders of magnitude between H = 2 and here.
 *
 * The tolerances below are still met with room to spare, because they are
 * relative to G(z) - G(0) and G(0) = V^2 / [2 pi A^2 (1 + 2H)] diverges as
 * H -> -1/2, growing the quantity being measured faster than the error grows.
 * That is the honest reading: the quadrature is working harder here, but the
 * observable is larger, so the error quoted on it is smaller.
 *
 * nu = H + 1/2 = 0.01 also makes the correlation function extremely cuspy - it
 * has already fallen to 5% of its full range by z = 1 - so there is no wide
 * transition region to sample whatever grid is chosen. */
static const gdab_case CASE_SHALLOW_DECAY = {
    .label = "H=-0.49 (slowest convergent tail)",
    .A = 10.0,
    .H = -0.49,
    .eta = 1e-4,
    .len = 12,
    .z = {1., 2., 5., 10., 25., 50., 100., 150., 200., 300., 400., 500.},
};

void setUp(void) {}
void tearDown(void) {}

/* G(0) for the GDAB model, written out rather than taken from the library so
 * that the reference the strategies are held to does not depend on the code
 * under test. test_analytical_gdab_tail_reaches_minus_G0 below checks that
 * this agrees with what compute_analytical_gdab itself converges to. */
static double gdab_G0(double A, double H, double eta) {
    double V = pow(2.0 * A, 3.0) * M_PI * sqrt(M_PI) * (tgamma(H + 1.5) / tgamma(H));
    return eta * eta * V * V / (2.0 * M_PI * A * A * (1.0 + 2.0 * H));
}

/* Runs one strategy over one parameter set and asserts every point is within
 * the strategy's relative tolerance of the analytical curve. */
static void check_strategy(const strategy_case *sc, const gdab_case *gc) {
    double analytical[MAX_LEN];
    double transform[MAX_LEN];
    double params[50] = {0};
    char message[256];

    params[0] = gc->A;
    params[1] = gc->H;
    params[2] = gc->eta;

    /* The form factor reads the same three parameters in the same order. */
    form_factor_ctx f_ctx = {.params = params};

    /* Cast away const on z: hankel_transform takes a mutable double *, but the
     * abscissae are only read. */
    double *z = (double *)gc->z;

    int status = compute_analytical_gdab(&params, z, analytical, gc->len);
    TEST_ASSERT_EQUAL_INT_MESSAGE(1, status, gc->label);

    /* Fixed_DE_Ogata needs f_max near the scale of the integrand's maximum;
     * the default 1.0 leaves it wrong by tens of percent on this form factor
     * while still reporting success. The other strategies ignore f_max. */
    strategy_params sp = {.n_eval = 250, .eps_rel = 1e-12, .f_max = 1e-4};

    status = hankel_transform(0, form_factor_g_dab, z, gc->len, (void *)&f_ctx, transform, sc->name,
                              sp);
    snprintf(message, sizeof(message), "%s on %s: transform failed", sc->name, gc->label);
    TEST_ASSERT_EQUAL_INT_MESSAGE(0, status, message);

    double G0 = gdab_G0(gc->A, gc->H, gc->eta);

    for (size_t i = 0; i < gc->len; i++) {
        double numerical = transform[i] / (2.0 * M_PI) - G0;

        snprintf(message, sizeof(message), "%s on %s at z=%g: rel err %.3e exceeds %.0e",
                 sc->name, gc->label, gc->z[i],
                 fabs(numerical - analytical[i]) / fabs(analytical[i]), sc->rel_tol);

        TEST_ASSERT_DOUBLE_WITHIN_MESSAGE(sc->rel_tol * fabs(analytical[i]), analytical[i],
                                          numerical, message);
    }
}

void test_strategies_match_analytical_gdab_integer_order(void) {
    for (size_t s = 0; s < N_STRATEGIES; s++) {
        check_strategy(&STRATEGIES[s], &CASE_INTEGER_ORDER);
    }
}

void test_strategies_match_analytical_gdab_real_order(void) {
    for (size_t s = 0; s < N_STRATEGIES; s++) {
        check_strategy(&STRATEGIES[s], &CASE_REAL_ORDER);
    }
}

void test_strategies_match_analytical_gdab_shallow_decay(void) {
    for (size_t s = 0; s < N_STRATEGIES; s++) {
        check_strategy(&STRATEGIES[s], &CASE_SHALLOW_DECAY);
    }
}

void test_analytical_gdab_tail_reaches_minus_G0(void) {
    /*
    G(z) - G(0) tends to -G(0) as the correlation dies out, so the constant
    subtracted from the numerical side above is pinned to the analytical curve
    itself and cannot drift from it unnoticed.
    */
    double params[50] = {10.0, 0.3, 1e-4};
    double z_far = 1e5;
    double tail;

    TEST_ASSERT_EQUAL_INT(1, compute_analytical_gdab(&params, &z_far, &tail, 1));
    TEST_ASSERT_DOUBLE_WITHIN(1e-12 * gdab_G0(10.0, 0.3, 1e-4), -gdab_G0(10.0, 0.3, 1e-4), tail);
}

void test_analytical_gdab_is_zero_at_z_zero(void) {
    /*
    G(0) - G(0) = 0. Worth pinning because the general expression is indefinite
    there: K_nu diverges while (u/2)^nu vanishes.
    */
    double params[50] = {10.0, 0.3, 1e-4};
    double z_zero = 0.0;
    double out = 1.0;

    TEST_ASSERT_EQUAL_INT(1, compute_analytical_gdab(&params, &z_zero, &out, 1));
    TEST_ASSERT_EQUAL_DOUBLE(0.0, out);
}

void test_analytical_gdab_rejects_out_of_range_parameters(void) {
    /*
    A = 0 has to be caught with the negatives: it divides by zero further down
    and would otherwise fill the output with NaN while reporting success.
    */
    double z = 15.0;
    double out = 0.0;
    char captured[256];

    double bad_A[50] = {0.0, 0.5, 1e-4};
    start_capture_stderr();
    int status = compute_analytical_gdab(&bad_A, &z, &out, 1);
    stop_capture_stderr(captured, sizeof(captured));
    TEST_ASSERT_EQUAL_INT_MESSAGE(-1, status, "A = 0 must be rejected");

    double negative_A[50] = {-10.0, 0.5, 1e-4};
    start_capture_stderr();
    status = compute_analytical_gdab(&negative_A, &z, &out, 1);
    stop_capture_stderr(captured, sizeof(captured));
    TEST_ASSERT_EQUAL_INT_MESSAGE(-1, status, "A < 0 must be rejected");

    double bad_H[50] = {10.0, -0.5, 1e-4};
    start_capture_stderr();
    status = compute_analytical_gdab(&bad_H, &z, &out, 1);
    stop_capture_stderr(captured, sizeof(captured));
    TEST_ASSERT_EQUAL_INT_MESSAGE(-1, status, "H = -1/2 must be rejected");

    double good[50] = {10.0, 0.5, 1e-4};
    double negative_z = -15.0;
    start_capture_stderr();
    status = compute_analytical_gdab(&good, &negative_z, &out, 1);
    stop_capture_stderr(captured, sizeof(captured));
    TEST_ASSERT_EQUAL_INT_MESSAGE(-1, status, "negative z must be rejected");
}

int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_strategies_match_analytical_gdab_integer_order);
    RUN_TEST(test_strategies_match_analytical_gdab_real_order);
    RUN_TEST(test_strategies_match_analytical_gdab_shallow_decay);
    RUN_TEST(test_analytical_gdab_tail_reaches_minus_G0);
    RUN_TEST(test_analytical_gdab_is_zero_at_z_zero);
    RUN_TEST(test_analytical_gdab_rejects_out_of_range_parameters);
    return UNITY_END();
}
