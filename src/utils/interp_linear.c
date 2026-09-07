// interp_linear.c
//
// Piecewise linear interpolation, the same API shape as interp_cubic so the
// two can be swapped.
//
// The arithmetic is written for three guarantees: exact at the nodes, never
// outside the bracketing pair of ordinates, and monotone in xi. The textbook
// one-liner delivers none of them reliably.

#include "interp_linear.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

struct linear_interp {
    double *x;
    double *y;
    size_t n;
};

/**
 * Blends f0 and f1 with weight t, which eval supplies in (0, 1) up to one
 * rounding of the quotient.
 *
 * Uses f0 + t*(f1 - f0) rather than the weighted average t*f1 + (1 - t)*f0.
 * Both are exact at t = 0, but only the correction form is monotone in t: the
 * weighted form's two products round independently, so a densely sampled curve
 * can wobble.
 *
 * Assumes f1 - f0 does not overflow, which needs ordinates near 1e308.
 */
static double blend(double f0, double f1, double t) {
    double y = f0 + t * (f1 - f0);

    /* The correction carries the sign of f1 - f0, so y can only ever drift
     * past the f1 end, including when t rounds up to 1 just inside a narrow
     * interval. This is the only thing holding the result inside [f0, f1]. */
    return (f1 > f0) ? fmin(y, f1) : fmax(y, f1);
}

/**
 * Index i with x[i] <= xi < x[i + 1]. The caller has already established
 * x[0] <= xi <= x[n - 1] and n >= 2, so such an i always exists.
 */
static size_t find_interval(const double *x, size_t n, double xi) {
    size_t lo = 0;
    size_t hi = n - 1;

    while (hi - lo > 1) {
        size_t mid = lo + (hi - lo) / 2;
        if (xi < x[mid]) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    return lo;
}

linear_interp_t *linear_interp_create(const double *x, const double *y, size_t n) {
    if (x == NULL || y == NULL || n < 2) {
        return NULL;
    }

    for (size_t j = 0; j < n; j++) {
        /* The > test is false for NaN, so it doubles as a NaN check; infinity
         * slips through it and needs the explicit isfinite. */
        if (j > 0 && !(x[j] > x[j - 1])) {
            return NULL;
        }
        if (!isfinite(x[j])) {
            return NULL;
        }
    }

    linear_interp_t *h = malloc(sizeof *h);
    if (h == NULL) {
        return NULL;
    }

    h->x = malloc(n * sizeof *h->x);
    h->y = malloc(n * sizeof *h->y);
    if (h->x == NULL || h->y == NULL) {
        free(h->x);
        free(h->y);
        free(h);
        return NULL;
    }

    memcpy(h->x, x, n * sizeof *h->x);
    memcpy(h->y, y, n * sizeof *h->y);
    h->n = n;
    return h;
}

double linear_interp_eval(const linear_interp_t *h, double xi) {
    if (h == NULL) {
        return NAN;
    }

    /* Negated so a NaN xi returns NaN rather than landing in some interval. */
    if (!(xi >= h->x[0] && xi <= h->x[h->n - 1])) {
        return NAN;
    }

    size_t i = find_interval(h->x, h->n, xi);

    /* The nodes are data, and these two tests are where the header's bit-exact
     * promise is kept -- the division below would only reproduce them to
     * within a rounding. Both ends are needed: find_interval binds xi to the
     * left of its interval, except at the top node, which lands on x[i + 1].
     */
    if (xi == h->x[i]) {
        return h->y[i];
    }
    if (xi == h->x[i + 1]) {
        return h->y[i + 1];
    }

    /* xi is strictly inside now, so t is in (0, 1) bar one rounding of the
     * quotient; a t that rounds up to 1 is absorbed by the clamp in blend. */
    double t = (xi - h->x[i]) / (h->x[i + 1] - h->x[i]);

    return blend(h->y[i], h->y[i + 1], t);
}

void linear_interp_destroy(linear_interp_t *h) {
    if (h == NULL) {
        return;
    }
    free(h->x);
    free(h->y);
    free(h);
}
