/**
 * @file interp_loglinear.h
 * @brief Log-space linear interpolation for positive tabulated data.
 *
 * Interpolates in log(y) space, returning exp(result), which gives exponential
 * (power-law) interpolation. Ideal for form factors and other data spanning
 * orders of magnitude: provides constant relative error across decades.
 *
 * Same API shape as @ref interp_linear.h and @ref interp_cubic.h so the three
 * are interchangeable wherever a table is interpolated.
 *
 * Compared with linear in real space: better relative accuracy over wide
 * ranges, but kinks at nodes remain (same as linear).
 *
 * @warning All ordinates must be strictly positive. @ref loglinear_interp_create
 *          rejects y <= 0. Do not pass @ref loglinear_interp_eval directly to
 *          @ref hankel_transform as a form factor; it returns NaN outside the
 *          tabulated range. See the
 *          <a href="../usage/tabulated_form_factors.html">Tabulated Form
 *          Factors</a> page.
 */

#ifndef INTERP_LOGLINEAR_H
#define INTERP_LOGLINEAR_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @typedef loglinear_interp_t
 * @brief Opaque handle to a built log-space interpolant.
 *
 * Only ever used through a pointer; do not dereference it or take its
 * `sizeof`. Same pattern as `FILE *` from `<stdio.h>`.
 */
typedef struct loglinear_interp loglinear_interp_t;

/**
 * @brief Builds a log-space linear interpolant through the given points.
 *
 * The input arrays are copied, so @p x and @p y may be freed as soon as this
 * returns.
 *
 * @param x  strictly increasing abscissae, at least @p n elements
 * @param y  positive ordinates, at least @p n elements
 * @param n  number of points; must be at least 2
 *
 * @return a handle to pass to @ref loglinear_interp_eval, or NULL if @p x or
 *         @p y is NULL, @p n is less than 2, @p x is not strictly increasing
 *         or holds a non-finite value, any @p y is non-positive or non-finite,
 *         or allocation failed. Every non-NULL return must eventually be passed
 *         to @ref loglinear_interp_destroy.
 */
loglinear_interp_t *loglinear_interp_create(const double *x, const double *y, size_t n);

/**
 * @brief Evaluates the interpolant at @p xi.
 *
 * Exact at the nodes: passing `x[i]` returns `y[i]` bit for bit. Between them
 * the result is guaranteed to be positive and monotone in @p xi (inherited
 * from monotone y data), so an interpolated intensity can never come back
 * negative or zero from positive data.
 *
 * @param h   handle from @ref loglinear_interp_create
 * @param xi  point at which to evaluate
 *
 * @return the interpolated value, or NaN if @p h is NULL, @p xi is NaN, or
 *         @p xi lies outside `[x[0], x[n-1]]`. This does not extrapolate.
 */
double loglinear_interp_eval(const loglinear_interp_t *h, double xi);

/**
 * @brief Frees an interpolant built by @ref loglinear_interp_create.
 *
 * @param h  handle to free; passing NULL is a safe no-op.
 */
void loglinear_interp_destroy(loglinear_interp_t *h);

#ifdef __cplusplus
}
#endif

#endif
