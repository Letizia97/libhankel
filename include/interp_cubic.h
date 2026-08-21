/**
 * @file interp_cubic.h
 * @brief Cubic interpolation, backed by Boost.Math's PCHIP interpolator.
 *
 * PCHIP is a cubic Hermite spline whose derivatives are chosen to preserve
 * monotonicity: it never overshoots the input data. That matters for
 * scattering curves, where a classical C2 spline can ring and produce
 * negative intensities. The price is that it is only C1 continuous and is
 * not exact for polynomials (a parabola comes back slightly low).
 *
 * Accuracy is worst in the first and last interval, where Boost falls back
 * to a one-sided O(h) estimate of the endpoint derivative.
 *
 * Typical use:
 * @code
 * cubic_interp_t *h = cubic_interp_create(x, y, n);
 * if (h == NULL) { ... }                 // n < 4, or x not increasing
 * for (i = 0; i < m; i++)
 *     out[i] = cubic_interp_eval(h, xi[i]);
 * cubic_interp_destroy(h);
 * @endcode
 *
 * @warning Do **not** pass @ref cubic_interp_eval directly to
 *          @ref hankel_transform as a form factor. The interpolator returns
 *          NaN outside the tabulated range, and every strategy evaluates the
 *          form factor far beyond it. See the
 *          <a href="../usage/tabulated_form_factors.html">Tabulated Form
 *          Factors</a> page, which explains what to do instead.
 */

/* Include guard: stops this file being pasted twice into one translation
 * unit. Not to be confused with the __cplusplus guards below, which do
 * something completely different. */

#ifndef INTERP_CUBIC_H
#define INTERP_CUBIC_H

#include <stddef.h>

/* Linkage guard. The implementation is C++ (interp_cubic.cpp), but this
 * header is included from C too. __cplusplus is only defined by a C++
 * compiler, so:
 *   - compiled as C++ -> extern "C" appears, symbols stay unmangled
 *   - compiled as C   -> it vanishes, which is required because
 *                        `extern "C"` is not valid C syntax
 * Without this, the C++ definition and the declaration seen here disagree
 * on linkage and the compiler rejects the file. */
#ifdef __cplusplus
extern "C" {
#endif

/* Opaque handle. `struct cubic_interp` is deliberately never defined in this
 * header -- its real definition lives in interp_cubic.cpp and holds a C++
 * Boost object that C could not parse. Callers may only hold and pass the
 * pointer, never dereference it or take sizeof.
 *
 * Because the type is incomplete, it can only ever be used through a pointer:
 * returning one by value is impossible, since the compiler cannot know its
 * size. Same pattern as FILE* from <stdio.h>. */

/**
 * @typedef cubic_interp_t
 * @brief Opaque handle to a built spline.
 *
 * The struct is deliberately incomplete: it holds a C++ object, so it can
 * only be used through a pointer. Do not dereference it or take its
 * `sizeof`. Same pattern as `FILE *` from `<stdio.h>`.
 */
typedef struct cubic_interp cubic_interp_t;

/**
 * @brief Builds a PCHIP spline through the given points.
 *
 * The input arrays are copied, so @p x and @p y may be freed or go out of
 * scope as soon as this returns.
 *
 * @param x  strictly increasing abscissae, at least @p n elements
 * @param y  ordinates, at least @p n elements
 * @param n  number of points; must be at least 4
 *
 * @return a handle to pass to @ref cubic_interp_eval, or NULL if @p x or
 *         @p y is NULL, @p n is less than 4, @p x is not strictly
 *         increasing, or allocation failed. Every non-NULL return must
 *         eventually be passed to @ref cubic_interp_destroy.
 */
cubic_interp_t *cubic_interp_create(const double *x, const double *y, size_t n);

/**
 * @brief Evaluates the spline at @p xi.
 *
 * `const` because evaluation does not modify the spline; this mirrors Boost,
 * whose `operator()` is a const member function.
 *
 * @param h   handle from @ref cubic_interp_create
 * @param xi  point at which to evaluate
 *
 * @return the interpolated value, or NaN if @p h is NULL or @p xi lies
 *         outside `[x[0], x[n-1]]`. Boost refuses to extrapolate.
 *
 * @warning The NaN-outside-range behaviour makes this function unsafe to use
 *          directly as a @ref form_factor_f. See the
 *          <a href="../usage/tabulated_form_factors.html">Tabulated Form
 *          Factors</a> page.
 */
double cubic_interp_eval(const cubic_interp_t *h, double xi);

/**
 * @brief Frees a spline built by @ref cubic_interp_create.
 *
 * @param h  handle to free; passing NULL is a safe no-op.
 */
void cubic_interp_destroy(cubic_interp_t *h);

#ifdef __cplusplus
}
#endif

#endif
