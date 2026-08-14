/* interp_cubic.h
 *
 * Cubic interpolation, backed by Boost.Math's PCHIP interpolator.
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
 * Usage:
 *     cubic_interp_t *h = cubic_interp_create(x, y, n);
 *     if (h == NULL) { ... }                 // n < 4, or x not increasing
 *     for (i = 0; i < m; i++)
 *         out[i] = cubic_interp_eval(h, xi[i]);
 *     cubic_interp_destroy(h);
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
typedef struct cubic_interp cubic_interp_t;

/* Builds the spline once. x must be strictly increasing and n >= 4.
 * Returns NULL on invalid input or allocation failure -- check it, the way
 * you would check fopen() or malloc(). */
cubic_interp_t *cubic_interp_create(const double *x, const double *y, size_t n);

/* Evaluates at xi. const because evaluation does not modify the spline;
 * this mirrors Boost, whose operator() is a const member function.
 *
 * Returns NaN if h is NULL or xi is outside [x[0], x[n-1]] -- Boost refuses
 * to extrapolate. Detect with isnan(). */
double cubic_interp_eval(const cubic_interp_t *h, double xi);

/* Frees the spline. Not const: this destroys the object rather than reading
 * it. Passing NULL is a safe no-op. */
void cubic_interp_destroy(cubic_interp_t *h);

#ifdef __cplusplus
}
#endif

#endif
