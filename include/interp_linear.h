/**
 * @file interp_linear.h
 * @brief Piecewise linear interpolation over a tabulated function.
 *
 * Deliberately the same shape as @ref interp_cubic.h -- create, eval, destroy,
 * an opaque handle, NaN outside the tabulated range -- so the two are
 * interchangeable wherever a table is interpolated.
 *
 * Compared with PCHIP: never overshoots at all rather than merely preserving
 * monotonicity, reproduces the tabulated values exactly, and is cheaper per
 * evaluation. The price is a kink at every node and O(h^2) error inside each
 * interval, against PCHIP's O(h^4).
 *
 * @warning Do **not** pass @ref linear_interp_eval directly to
 *          @ref hankel_transform as a form factor. Like the cubic
 *          interpolator it returns NaN outside the tabulated range, and every
 *          strategy evaluates the form factor far beyond it. See the
 *          <a href="../usage/tabulated_form_factors.html">Tabulated Form
 *          Factors</a> page.
 */

#ifndef INTERP_LINEAR_H
#define INTERP_LINEAR_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @typedef linear_interp_t
 * @brief Opaque handle to a built interpolant.
 *
 * Only ever used through a pointer; do not dereference it or take its
 * `sizeof`. Same pattern as `FILE *` from `<stdio.h>`.
 */
typedef struct linear_interp linear_interp_t;

/**
 * @brief Builds a piecewise linear interpolant through the given points.
 *
 * The input arrays are copied, so @p x and @p y may be freed as soon as this
 * returns.
 *
 * @param x  strictly increasing abscissae, at least @p n elements
 * @param y  ordinates, at least @p n elements
 * @param n  number of points; must be at least 2
 *
 * @return a handle to pass to @ref linear_interp_eval, or NULL if @p x or
 *         @p y is NULL, @p n is less than 2, @p x is not strictly increasing
 *         or holds a non-finite value, or allocation failed. Every non-NULL
 *         return must eventually be passed to @ref linear_interp_destroy.
 */
linear_interp_t *linear_interp_create(const double *x, const double *y, size_t n);

/**
 * @brief Evaluates the interpolant at @p xi.
 *
 * Exact at the nodes: passing `x[i]` returns `y[i]` bit for bit. Between them
 * the result is guaranteed to lie within `[y[i], y[i+1]]` and to be monotone
 * in @p xi, so an interpolated intensity can never come back negative from
 * non-negative data.
 *
 * @param h   handle from @ref linear_interp_create
 * @param xi  point at which to evaluate
 *
 * @return the interpolated value, or NaN if @p h is NULL, @p xi is NaN, or
 *         @p xi lies outside `[x[0], x[n-1]]`. This does not extrapolate.
 */
double linear_interp_eval(const linear_interp_t *h, double xi);

/**
 * @brief Frees an interpolant built by @ref linear_interp_create.
 *
 * @param h  handle to free; passing NULL is a safe no-op.
 */
void linear_interp_destroy(linear_interp_t *h);

#ifdef __cplusplus
}
#endif

#endif
