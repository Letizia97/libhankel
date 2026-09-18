// interp_loglinear.c
//
// Log-space linear interpolation. Wraps interp_linear.c, interpolating in
// log(y) space and returning exp(result). Gives exponential (power-law)
// interpolation suitable for form factors spanning orders of magnitude.

#include "interp_loglinear.h"

#include <math.h>
#include <stdlib.h>

#include "interp_linear.h"

struct loglinear_interp {
    linear_interp_t *lin;
};

loglinear_interp_t *loglinear_interp_create(const double *x, const double *y, size_t n) {
    if (x == NULL || y == NULL || n < 2) {
        return NULL;
    }

    // Validate x: strictly increasing and finite.
    for (size_t j = 0; j < n; j++) {
        if (j > 0 && !(x[j] > x[j - 1])) {
            return NULL;
        }
        if (!isfinite(x[j])) {
            return NULL;
        }
    }

    // Validate y: all positive and finite.
    for (size_t j = 0; j < n; j++) {
        if (y[j] <= 0.0 || !isfinite(y[j])) {
            return NULL;
        }
    }

    // Allocate the loglinear handle.
    loglinear_interp_t *h = malloc(sizeof *h);
    if (h == NULL) {
        return NULL;
    }

    // Convert y to log(y) in a temporary buffer and build linear interp.
    // same as malloc(n * sizeof(double)) as dereferencing log_y gives a double
    double *log_y = malloc(n * sizeof *log_y);
    if (log_y == NULL) {
        free(h);
        return NULL;
    }

    for (size_t j = 0; j < n; j++) {
        log_y[j] = log(y[j]);
    }

    h->lin = linear_interp_create(x, log_y, n);
    free(log_y);

    if (h->lin == NULL) {
        free(h);
        return NULL;
    }

    return h;
}

double loglinear_interp_eval(const loglinear_interp_t *h, double xi) {
    if (h == NULL) {
        return NAN;
    }

    double log_result = linear_interp_eval(h->lin, xi);

    if (!isfinite(log_result)) {
        return NAN;
    }

    return exp(log_result);
}

void loglinear_interp_destroy(loglinear_interp_t *h) {
    if (h == NULL) {
        return;
    }
    linear_interp_destroy(h->lin);
    free(h);
}
