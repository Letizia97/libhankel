// interp_linear.c
//
// Piecewise linear interpolation, the same API shape as interp_cubic so the
// two can be swapped. Plain C: unlike PCHIP there is nothing here worth
// pulling Boost in for.
//
// The arithmetic is written for three guarantees, in this order: exact at the
// nodes, never outside the bracketing pair of ordinates, and monotone in xi.
// The textbook one-liner delivers none of them reliably.

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
 * Blends f0 and f1 with weight t in [0, 1], the numerically careful way.
 *
 * The two algebraically identical forms fail in opposite places. Given
 * f0 + t*(f1 - f0), the difference can overflow and cancel when the ordinates
 * straddle zero, and t = 1 does not return f1 exactly. Given
 * (1 - t)*f0 + t*f1, both endpoints are exact but the result is not monotone
 * in t, so a curve sampled densely can wobble.
 *
 * So: take the weighted form only where it is needed, pin the endpoint, and
 * clamp what remains. Follows the reasoning behind C++20's std::lerp, which
 * this project cannot use - it is C++17, and this file is C.
 */
static double blend(double f0, double f1, double t) {
    /* Straddling zero is the only case where f1 - f0 can overflow, and the
     * worst case for cancellation. Here both terms lie between f0 and f1, so
     * the weighted form cannot leave the interval and is exact at both ends. */
    if ((f0 <= 0.0 && f1 >= 0.0) || (f0 >= 0.0 && f1 <= 0.0)) {
        return t * f1 + (1.0 - t) * f0;
    }

    /* Same sign, so f0 + t*(f1 - f0) is the accurate form: it perturbs f0 by a
     * correction rather than differencing two comparable products. It is exact
     * at t = 0 for free, but t = 1 has to be pinned by hand. */
    if (t == 1.0) {
        return f1;
    }

    double y = f0 + t * (f1 - f0);

    /* Belt and braces. No case where this fires turned up in 10^8 same-sign
     * pairs, but the bound is not proven and it costs one fmin. Only the f1
     * end can overshoot, since the correction never changes sign. */
    return (f1 > f0) ? fmin(y, f1) : fmax(y, f1);
}

/**
 * Index i with x[i] <= xi < x[i + 1]. The caller has already established
 * x[0] <= xi <= x[n - 1] and n >= 2, so such an i always exists.
 */
static size_t find_interval(const double *x, size_t n, double xi) {
    size_t lo = 0;
    size_t hi = n - 1;

    /* Invariant x[lo] <= xi < x[hi], narrowed until the two are adjacent. */
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
        /* The strictly-increasing test is false for NaN, so it doubles as a
         * NaN check on x. It misses infinity, which would make every width
         * infinite and every t zero, hence the explicit isfinite. */
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
        /* free(NULL) is a no-op, so this covers either allocation failing. */
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

    /* Negated so a NaN xi fails the test and returns NaN rather than being
     * silently placed in an interval. */
    if (!(xi >= h->x[0] && xi <= h->x[h->n - 1])) {
        return NAN;
    }

    size_t i = find_interval(h->x, h->n, xi);

    /* Answered before any arithmetic: the nodes are data, and the division
     * below would only reproduce them to within a rounding error. */
    if (xi == h->x[i]) {
        return h->y[i];
    }
    if (xi == h->x[i + 1]) {
        return h->y[i + 1];
    }

    /* Both differences are exact whenever the nodes are within a factor of two
     * of each other (Sterbenz), which is the common case for a log-spaced
     * table. The quotient is the only rounding on this path. */
    double t = (xi - h->x[i]) / (h->x[i + 1] - h->x[i]);

    /* That single rounding is enough to give t = 1 + eps on a narrow interval,
     * which would extrapolate a hair past the node. */
    if (t < 0.0) {
        t = 0.0;
    } else if (t > 1.0) {
        t = 1.0;
    }

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
