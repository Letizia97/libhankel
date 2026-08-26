// tabulated_ff.c
//
// A form_factor_f backed by a table of points rather than a formula.
//
// interp_cubic interpolates; this adds the tails. The strategies sample q over
// tens of decades, so every call falls outside the table, where
// cubic_interp_eval returns NaN. A tail that would not converge is refused at
// build time.

#include "tabulated_ff.h"

#include "interp_cubic.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Endpoints are cached because most calls need them. */
struct tabulated_ff {
    cubic_interp_t *spline;
    double q_lo; /* first tabulated q                                  */
    double q_hi; /* last tabulated q                                   */
    double f_lo; /* form factor at q_lo, held flat below it            */
    double f_hi; /* form factor at q_hi, where the high-q tail joins   */
    tabulated_tail_t tail;
    double exponent; /* the p of the power-law tail; unused for _ZERO  */
};

/* Points from the high-q end used to fit the slope: a quarter of the table,
 * clamped. Below four is not a fit; above forty reaches into the pre-asymptotic
 * curve and biases p towards zero, the direction that diverges. */
#define FIT_MIN_POINTS 4
#define FIT_MAX_POINTS 40

/* Below this the transform does not converge: the integrand goes like
 * q^(1/2 - p), so the integral only settles for p > 3/2. */
#define MIN_CONVERGENT_EXPONENT 1.5

/**
 * Rejects an unusable table before anything is allocated. Returns -13.
 *
 * The strictly-increasing test doubles as a NaN check on q, since any comparison
 * with NaN is false. It misses infinity, hence the explicit isfinite calls.
 */
static int validate_table(const double *q, const double *f, size_t n) {
    if (q == NULL || f == NULL) {
        fprintf(stderr, "Error: tabulated form factor needs non-NULL q and f\n");
        return -13;
    }
    if (n < FIT_MIN_POINTS) {
        fprintf(stderr, "Error: tabulated form factor needs at least %d points, got %zu\n",
                FIT_MIN_POINTS, n);
        return -13;
    }
    if (!(q[0] > 0)) {
        fprintf(stderr, "Error: tabulated q must be greater than zero, got q[0] = %g\n", q[0]);
        return -13;
    }
    for (size_t j = 0; j < n; j++) {
        if (j > 0 && !(q[j] > q[j - 1])) {
            fprintf(stderr, "Error: tabulated q must be strictly increasing, but q[%zu] = %g is "
                            "not greater than q[%zu] = %g\n",
                    j, q[j], j - 1, q[j - 1]);
            return -13;
        }
        /* An infinite q_hi never satisfies `q >= q_hi`, so the tail is
         * unreachable and every eval returns NaN under status 0. */
        if (!isfinite(q[j])) {
            fprintf(stderr, "Error: tabulated q must be finite, but q[%zu] = %g\n", j, q[j]);
            return -13;
        }
        /* A non-finite f builds a spline happily, then poisons every eval. */
        if (!isfinite(f[j])) {
            fprintf(stderr, "Error: tabulated f must be finite, but f[%zu] = %g\n", j, f[j]);
            return -13;
        }
    }
    return 0;
}

/**
 * Fits the high-q slope by least squares of log f against log q, since
 * f = A q^-p is a line of gradient -p in log-log. Only the high-q end is
 * usable: the low-q plateau has a near-zero slope, which diverges.
 *
 * Points with f <= 0 are skipped rather than rejected, because
 * background-subtracted data goes negative exactly where this fit looks.
 *
 * Writes p to *out_p and returns 0, or -14.
 */
static int fit_exponent(const double *q, const double *f, size_t n, double *out_p) {
    size_t count = n / 4;
    if (count < FIT_MIN_POINTS) {
        count = FIT_MIN_POINTS;
    }
    if (count > FIT_MAX_POINTS) {
        count = FIT_MAX_POINTS;
    }
    if (count > n) {
        count = n;
    }

    double sum_x = 0.0, sum_y = 0.0, sum_xx = 0.0, sum_xy = 0.0;
    size_t used = 0;

    for (size_t j = n - count; j < n; j++) {
        if (!(f[j] > 0)) {
            continue;
        }
        double x = log(q[j]);
        double y = log(f[j]);
        sum_x += x;
        sum_y += y;
        sum_xx += x * x;
        sum_xy += x * y;
        used++;
    }

    /* Two points fix a gradient with no redundancy, so one noisy value would
     * set the tail for the whole transform. */
    if (used < 3) {
        fprintf(stderr,
                "Error: cannot fit a high-q exponent, only %zu of the last %zu tabulated points "
                "are positive. Pass an explicit exponent instead.\n",
                used, count);
        return -14;
    }

    double m = (double)used;
    double denominator = m * sum_xx - sum_x * sum_x;
    double gradient = (m * sum_xy - sum_x * sum_y) / denominator;

    /* Checking the result rather than the denominator catches a q range too
     * narrow to fit, and any overflow, at once. */
    if (!isfinite(gradient)) {
        fprintf(stderr, "Error: high-q exponent fit did not produce a finite gradient\n");
        return -14;
    }

    /* f ~ q^gradient, and the tail is written as q^-p. */
    *out_p = -gradient;
    return 0;
}

int tabulated_ff_create(const double *q, const double *f, size_t n, tabulated_tail_t tail,
                        double exponent, tabulated_ff_t **out) {
    if (out == NULL) {
        fprintf(stderr, "Error: tabulated_ff_create needs a non-NULL out pointer\n");
        return -13;
    }
    /* Cleared up front so every failure below can simply return. */
    *out = NULL;

    int status = validate_table(q, f, n);
    if (status != 0) {
        return status;
    }

    if (tail != TABULATED_TAIL_POWER_LAW && tail != TABULATED_TAIL_ZERO) {
        fprintf(stderr, "Error: unknown tail policy %d\n", (int)tail);
        return -13;
    }

    /* Settled before allocating; the spline costs tens of microseconds. */
    double p = 0.0;
    if (tail == TABULATED_TAIL_POWER_LAW) {
        if (exponent == 0.0) {
            /* Zero is free to mean "fit it": p = 0 is the flat tail, which
             * provably does not converge. */
            status = fit_exponent(q, f, n, &p);
            if (status != 0) {
                return status;
            }
        } else {
            p = exponent;
        }

        if (!isfinite(p) || p <= MIN_CONVERGENT_EXPONENT) {
            fprintf(stderr,
                    "Error: high-q exponent %g does not converge, it must be greater than %g. If "
                    "this was fitted, the table most likely stops before the asymptotic regime.\n",
                    p, MIN_CONVERGENT_EXPONENT);
            return -14;
        }
    }

    tabulated_ff_t *t = malloc(sizeof *t);
    if (t == NULL) {
        fprintf(stderr, "Error: failed to allocate tabulated form factor\n");
        return -3;
    }

    /* Copies q and f, which is what frees the caller from keeping them alive. */
    t->spline = cubic_interp_create(q, f, n);
    if (t->spline == NULL) {
        /* The table is already validated, so this is an allocation failure. */
        fprintf(stderr, "Error: failed to build the interpolating spline\n");
        free(t);
        return -3;
    }

    t->q_lo = q[0];
    t->q_hi = q[n - 1];
    t->f_lo = f[0];
    t->f_hi = f[n - 1];
    t->tail = tail;
    t->exponent = p;

    *out = t;
    return 0;
}

double tabulated_ff_eval(double q, void *ctx) {
    const tabulated_ff_t *t = (const tabulated_ff_t *)ctx;

    if (t == NULL) {
        return NAN;
    }

    /* Inclusive comparisons, so the endpoints come back exactly. A NaN q fails
     * both and falls through to the spline, which returns NaN - the right
     * answer. */
    if (q <= t->q_lo) {
        return t->f_lo;
    }

    if (q >= t->q_hi) {
        /* The last point is data, not tail: a truncated tail would otherwise
         * discard a measured value here. */
        if (q == t->q_hi) {
            return t->f_hi;
        }
        if (t->tail == TABULATED_TAIL_ZERO) {
            return 0.0;
        }
        /* A ratio rather than A * pow(q, -p), so the tail meets the spline
         * exactly at q_hi whatever p turned out to be. */
        return t->f_hi * pow(q / t->q_hi, -t->exponent);
    }

    return cubic_interp_eval(t->spline, q);
}

void tabulated_ff_destroy(tabulated_ff_t *t) {
    if (t == NULL) {
        return;
    }
    /* Never NULL: create frees t rather than hand back a handle without one. */
    cubic_interp_destroy(t->spline);
    free(t);
}
