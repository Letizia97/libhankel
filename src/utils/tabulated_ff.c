// tabulated_ff.c
//
// A form_factor_f backed by a table of points rather than a formula.
//
// The interpolation itself is done by interp_cubic; everything here exists
// because a spline alone is not enough. cubic_interp_eval returns NaN outside
// the tabulated range, and the strategies evaluate the form factor over tens
// of decades in q, so a bare spline is asked for values it cannot supply on
// every single call. This wraps it with the two tail rules that make it
// defined everywhere, and refuses to build one whose tail would not converge.

#include "tabulated_ff.h"

#include "interp_cubic.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* The handle the header keeps opaque. Both endpoint values are cached rather
 * than read back off the spline: they are needed on the majority of calls,
 * and asking the spline for its own boundary would be both slower and, at
 * exactly x[0] and x[n-1], the one place its range check could bite. */
struct tabulated_ff {
    cubic_interp_t *spline;
    double q_lo; /* first tabulated q                                  */
    double q_hi; /* last tabulated q                                   */
    double f_lo; /* form factor at q_lo, held flat below it            */
    double f_hi; /* form factor at q_hi, where the high-q tail joins   */
    tabulated_tail_t tail;
    double exponent; /* the p of the power-law tail; unused for _ZERO  */
};

/* Number of points from the high-q end of the table - that is, the last
 * entries of q and f - used to fit the high-q slope.
 *
 * A quarter of the table, so a short one still has something to work with,
 * bounded at both ends: fewer than four points is not a fit, and more than
 * forty starts reaching down into the part of the curve that has not yet
 * reached its asymptote, which biases the slope towards zero - the dangerous
 * direction, since that is the one that fails to converge. */
#define FIT_MIN_POINTS 4
#define FIT_MAX_POINTS 40

/* Below this the transform does not converge: the integrand goes like
 * q^(1/2 - p), so the integral only settles for p > 3/2. */
#define MIN_CONVERGENT_EXPONENT 1.5

/**
 * Rejects a table that cannot be used, before anything has been allocated.
 *
 * Returns -13 on any problem, matching tabulated_ff_create's contract.
 *
 * Note that the strictly-increasing test doubles as a NaN check on q: a
 * comparison against NaN is false, so `!(q[j] > q[j - 1])` catches it. The
 * same trick does not work on f, which is only ever read and never compared,
 * so that one needs an explicit isfinite.
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
        /* A NaN here would build a spline quite happily and then poison every
         * evaluation that touched it, which is exactly the silent failure
         * this whole type exists to prevent. */
        if (!isfinite(f[j])) {
            fprintf(stderr, "Error: tabulated f must be finite, but f[%zu] = %g\n", j, f[j]);
            return -13;
        }
    }
    return 0;
}

/**
 * Fits the high-q slope from the last points in the table.
 *
 * Only the high-q end can be used. A power law is an asymptotic statement, so
 * it holds where q is largest; the other end of the table is the low-q
 * plateau, whose slope is close to zero and therefore divergent.
 *
 * A power law f = A q^-p is a straight line of gradient -p once both axes are
 * logged, so this is an ordinary least-squares fit of log f against log q.
 *
 * Points with f <= 0 are skipped rather than treated as an error. Real data
 * that has had a background subtracted goes negative at high q, precisely
 * where this fit looks, and log of that is NaN - which would not fail loudly,
 * it would quietly poison every sum below and produce a NaN gradient.
 *
 * Writes the fitted p to *out_p and returns 0, or returns -14 if too few
 * usable points survived or the fit did not produce a finite gradient.
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

    /* Two points would define a gradient, but with no redundancy at all a
     * single noisy value sets the tail for the entire transform. Three is the
     * smallest fit that can be wrong in a visible way. */
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

    /* The abscissae are distinct, so the denominator is only at risk from
     * rounding on a very short, very narrow q range. Checking the result
     * rather than the denominator catches that and any overflow at once. */
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
    /* Cleared up front so that every failure below can simply return, and the
     * caller is never left holding a stale pointer it might try to destroy. */
    *out = NULL;

    int status = validate_table(q, f, n);
    if (status != 0) {
        return status;
    }

    if (tail != TABULATED_TAIL_POWER_LAW && tail != TABULATED_TAIL_ZERO) {
        fprintf(stderr, "Error: unknown tail policy %d\n", (int)tail);
        return -13;
    }

    /* Settled before anything is allocated: a table whose tail diverges is
     * not going to become usable later, and building a spline for it would be
     * tens of microseconds spent on a call that cannot succeed. */
    double p = 0.0;
    if (tail == TABULATED_TAIL_POWER_LAW) {
        if (exponent == 0.0) {
            /* Zero is free to mean "fit it" because it is never a legitimate
             * value: p = 0 is the flat tail, which is the one option that
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

    /* Copies q and f, which is what frees the caller from keeping them
     * alive - see the header. */
    t->spline = cubic_interp_create(q, f, n);
    if (t->spline == NULL) {
        /* The table has already been validated, so this is not a data
         * problem: the allocation inside the spline failed. */
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

    /* Both comparisons are inclusive, so the two endpoints are answered from
     * the cached values and never reach the spline. That is deliberate: it
     * returns the tabulated value exactly, and it keeps the boundary of the
     * spline's own range check out of the hot path entirely.
     *
     * A NaN q fails both tests and falls through to the spline, whose range
     * check rejects it and returns NaN - which is the right answer, and the
     * only input for which this function is not total. */
    if (q <= t->q_lo) {
        return t->f_lo;
    }

    if (q >= t->q_hi) {
        if (t->tail == TABULATED_TAIL_ZERO) {
            return 0.0;
        }
        /* Written as a ratio rather than as A * pow(q, -p) so the tail meets
         * the spline exactly at q_hi, where the ratio is 1. No step at the
         * join, whatever p turned out to be. */
        return t->f_hi * pow(q / t->q_hi, -t->exponent);
    }

    return cubic_interp_eval(t->spline, q);
}

void tabulated_ff_destroy(tabulated_ff_t *t) {
    if (t == NULL) {
        return;
    }
    /* cubic_interp_destroy is itself NULL-safe, but t->spline cannot be NULL
     * here anyway: create frees t and returns rather than handing back a
     * handle with no spline in it. */
    cubic_interp_destroy(t->spline);
    free(t);
}
