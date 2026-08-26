#ifndef TABULATED_FF_H
#define TABULATED_FF_H

#include <stddef.h>

/*
Usable from C++ as well as C: without this guard a C++ compiler would mangle
the three function names and the link against the C-compiled object would
fail. The mirror image of the guard in interp_cubic.h, where the
implementation is C++ and the caller C.
*/
#ifdef __cplusplus
extern "C" {
#endif

/**
 * @enum tabulated_tail_t
 * @brief How the form factor behaves above the last tabulated point.
 *
 * Not an edge case: the strategies sample @f$ q @f$ over tens of decades, so
 * this is hit on every call. See the
 * <a href="../usage/tabulated_form_factors.html">Tabulated form factors</a>
 * page for the measured ranges.
 *
 * Only the high-@f$ q @f$ side is a choice. Below the table the value is held
 * flat at the first tabulated point - the physical @f$ q \to 0 @f$ plateau,
 * harmless because the integrand's factor of @f$ q @f$ suppresses that region.
 */
typedef enum {
    /**
     * @brief Continue as @f$ f(q) = f(q_{n-1})\,(q/q_{n-1})^{-p} @f$.
     *
     * The integrand then goes like @f$ q^{1/2-p} @f$, so the transform
     * converges only for @f$ p > 3/2 @f$; @ref tabulated_ff_create rejects
     * anything at or below that. Porod's law gives @f$ p = 4 @f$ for a sharp
     * interface. Continuous at the join.
     */
    TABULATED_TAIL_POWER_LAW = 0,

    /**
     * @brief Drop to zero above the last tabulated point.
     *
     * Always converges, but the discontinuity puts ringing into the result.
     * Correct when the data really is zero beyond the table - an instrumental
     * @f$ q @f$ cutoff, say - and wrong when the table simply ran out.
     */
    TABULATED_TAIL_ZERO = 1
} tabulated_tail_t;

/**
 * @typedef tabulated_ff_t
 * @brief Opaque handle to a form factor defined by a table of points.
 *
 * A PCHIP spline through the data plus the tail rules that make it defined for
 * every @f$ q > 0 @f$. Hidden so the tail bookkeeping cannot be tampered with.
 */
typedef struct tabulated_ff tabulated_ff_t;

/**
 * @brief Builds a form factor from tabulated data.
 *
 * @p q and @p f are copied, so the caller may free or reuse them as soon as
 * this returns.
 *
 * @param q         strictly increasing, strictly positive abscissae
 * @param f         form factor value at each entry of @p q
 * @param n         number of points; at least 4, as a cubic through fewer is
 *                  underdetermined
 * @param tail      high-@f$ q @f$ behaviour, see @ref tabulated_tail_t
 * @param exponent  the @f$ p @f$ of @ref TABULATED_TAIL_POWER_LAW. Pass 0 to
 *                  fit it from the last points of the table by least squares
 *                  in log-log space. Ignored for @ref TABULATED_TAIL_ZERO.
 * @param out       receives the new handle on success, NULL on every failure
 *
 * @return 0 on success, or
 *         -3  if allocation failed,
 *         -13 if @p q, @p f or @p out is NULL, @p n is below 4, @p q is not
 *             strictly increasing and positive, @p q or @p f contains a value
 *             that is not finite, or @p tail is not one of
 *             @ref tabulated_tail_t, or
 *         -14 if the power-law tail is unusable: the exponent is at or below
 *             3/2, supplied or fitted, so the transform would not converge, or
 *             too few of the points it would have been fitted from are
 *             positive.
 *
 * @warning A -14 from a fitted exponent usually means the table stops before
 *          the asymptotic regime, not that the physics diverges. The fit can
 *          only see the slope it was given: a table truncated at
 *          @f$ q_{max}\xi = 1 @f$ yields @f$ p \approx 1.2 @f$ even when the
 *          true asymptote is 4. Extend the table, or pass the exponent you
 *          know applies.
 */
int tabulated_ff_create(const double *q, const double *f, size_t n, tabulated_tail_t tail,
                        double exponent, tabulated_ff_t **out);

/**
 * @brief Evaluates the tabulated form factor.
 *
 * Shaped as a @ref form_factor_f, so a handle from @ref tabulated_ff_create
 * goes straight to @ref hankel_transform as its @c f_ctx argument, with no
 * adapter to write.
 *
 * Unlike @ref cubic_interp_eval this is total: every @p q gets a finite
 * number, from the spline inside the table and the tail rules outside it. Both
 * tabulated endpoints are returned exactly, whichever tail is in force, since
 * the last point is data - a truncated tail begins above it, not at it.
 *
 * @param q    value at which to evaluate. Values at or below zero give the
 *             first tabulated value, continuing the low-@f$ q @f$ plateau.
 * @param ctx  handle from @ref tabulated_ff_create
 *
 * @return the form factor at @p q, or NaN if @p ctx is NULL
 */
double tabulated_ff_eval(double q, void *ctx);

/**
 * @brief Releases a handle from @ref tabulated_ff_create.
 *
 * Does nothing if @p t is NULL, so cleanup paths can call it unconditionally,
 * exactly as free() allows.
 *
 * @param t  handle to destroy
 */
void tabulated_ff_destroy(tabulated_ff_t *t);

#ifdef __cplusplus
}
#endif

#endif // TABULATED_FF_H
