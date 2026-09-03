#include "src/utils/analytical_form_factors.h"

#include <math.h>
#include <stddef.h>
#include <stdio.h>

#include "libhankel.h"

#include "../src/utils/sf_functions.h"
#include "../src/utils/boost_bessel_wrapper.h"

/**
 * @brief Analytical SESANS correlation function of a homogeneous sphere.
 *
 * Closed form of @f$ \eta^2 [G(z) - G(0)] @f$ for a sphere of radius @f$ R @f$,
 * that is, the 2D projection
 * @f$ \frac{1}{2\pi}\int_0^\infty I(q) J_0(qz)\, q \, dq @f$ of
 * @ref form_factor_sphere. Same role as @ref compute_analytical_gdab: a
 * reference against which the numerical Hankel strategies are checked.
 *
 * With @f$ \xi = z / 2R @f$,
 * @f[
 *   G(z) = \pi \eta^2 R^4 \left[ \sqrt{1-\xi^2}\,(2+\xi^2)
 *          + \xi^2 (4-\xi^2) \ln\frac{\xi}{1+\sqrt{1-\xi^2}} - 2 \right] ,
 * @f]
 * which is 0 at @f$ z = 0 @f$ and decreases to @f$ -2\pi \eta^2 R^4 @f$ at
 * @f$ z = 2R @f$. Two points in the sphere cannot be further apart than a
 * diameter, so the curve is exactly flat above @f$ 2R @f$.
 *
 * The logarithm is written with @f$ \xi @f$ in the numerator so that its
 * argument stays bounded as @f$ \xi \to 1 @f$; the divergence at
 * @f$ \xi \to 0 @f$ is harmless because the @f$ \xi^2 @f$ factor kills it, and
 * a small floor keeps `log(0)` out of the expression at @f$ z = 0 @f$ exactly.
 *
 * @param params  pointer to the parameter array: @p (*params)[0] is the radius
 *                @f$ R @f$, @p (*params)[1] the scattering length density
 *                contrast @f$ \eta @f$. Later entries are unused.
 * @param arr_z   spin-echo lengths, at least @p n elements
 * @param G       receives @f$ \eta^2 [G(z) - G(0)] @f$, at least @p n elements
 * @param n       number of points
 *
 * @return 1 if at least one input fell inside @f$ 0 \le z < 2R @f$, 0
 *         otherwise, in which case every element of @p G is left at zero.
 *
 * @note Points at @f$ z \ge 2R @f$ are filled with the last value computed
 *       inside the domain rather than with the exact plateau
 *       @f$ -2\pi \eta^2 R^4 @f$, so their accuracy depends on how finely
 *       @p arr_z samples the approach to @f$ 2R @f$: a grid whose last
 *       in-domain point is @f$ z = 19 @f$ with @f$ R = 10 @f$ puts the plateau
 *       0.16% low. Assigning @f$ -2\pi \eta^2 R^4 @f$ directly would remove
 *       the grid dependence.
 */
double compute_analytical_spheres(double (*params)[50], const double *arr_z, double *G, size_t n) {

    const double R = (*params)[0];
    const double ETA = (*params)[1];

    const double eps = 1e-15;
    const double pref = ETA * ETA * M_PI * pow(R, 4.0);

    double last_G_val = 0.0;
    int have_last = 0;

    for (size_t i = 0; i < n; i++) {
        double z = arr_z[i];
        double xi = z / (2.0 * R);

        // Default
        G[i] = 0.0;

        // Valid domain: 0 <= xi < 1
        if (xi >= 0.0 && xi < 1.0) {
            // printf("Entered block A\n");

            double xi2 = xi * xi;
            double sqrt_term = sqrt(1.0 - xi2);

            // First term: sqrt(1 - xi^2) * (2 + xi^2)
            double term1 = sqrt_term * (2.0 + xi2);

            // Log term
            double log_arg = xi / (1.0 + sqrt_term);
            if (log_arg < eps) {
                log_arg = eps;
            }

            double term_log = log(log_arg);
            double term2 = xi2 * (4.0 - xi2) * term_log;

            double G_val = pref * (term1 + term2 - 2.0);

            G[i] = G_val;
            last_G_val = G_val;
            have_last = 1;
        }
    }

    // Ensure G does not suddenly go to zero once z > 2R
    if (have_last) {
        for (size_t i = 0; i < n; i++) {
            double xi = arr_z[i] / (2.0 * R);
            if (!(xi >= 0.0 && xi < 1.0)) {
                G[i] = last_G_val;
            }
        }
        return 1;
    }

    return 0;
}

/**
 * @brief Analytical SESANS correlation function of the GDAB model.
 *
 * Closed form of the quantity a SESANS measurement is sensitive to,
 * @f$ \eta^2 [G(z) - G(0)] @f$, for the Generalized Debye-Anderson-Brumberger
 * (Whittle-Matern) model. It is the exact 2D projection
 * @f$ G(z) = \frac{1}{2\pi}\int_0^\infty I(q) J_0(qz)\, q \, dq @f$
 * of @ref form_factor_g_dab, so it is the reference against which the
 * numerical Hankel strategies are checked.
 *
 * With @f$ u = z/A @f$ and @f$ \nu = H + 1/2 @f$,
 * @f[
 *   G(z) = \frac{V^2}{2\pi A^2 \Gamma(\nu+1)} \left(\frac{u}{2}\right)^{\nu}
 *          K_{\nu}(u), \qquad
 *   V = (2A)^3 \pi^{3/2} \frac{\Gamma(H+3/2)}{\Gamma(H)} ,
 * @f]
 * whose @f$ z \to 0 @f$ limit is @f$ G(0) = V^2 / [2\pi A^2 (1+2H)] @f$.
 * The result is therefore zero at @f$ z = 0 @f$, negative elsewhere, and
 * tends to @f$ -\eta^2 G(0) @f$ as @f$ z \to \infty @f$.
 *
 * @param params  pointer to the parameter array: @p (*params)[0] is the
 *                correlation length @f$ A > 0 @f$, @p (*params)[1] the Hurst
 *                exponent @f$ H > -1/2 @f$, @p (*params)[2] the scattering
 *                length density contrast @f$ \eta @f$. Later entries are
 *                unused.
 * @param arr_z   spin-echo lengths, at least @p n elements, all non-negative
 * @param out     receives @f$ \eta^2 [G(z) - G(0)] @f$, at least @p n elements
 * @param n       number of points
 *
 * @return 1 on success, -1 if @p A, @p H or any element of @p arr_z is out of
 *         range. On failure @p out may already be partially written.
 *
 * @note @f$ \Gamma(H) @f$ in @f$ V @f$ diverges at @f$ H = 0 @f$, so the model
 *       amplitude vanishes there. That normalization is inherited from
 *       @ref form_factor_g_dab and is not a defect of this routine.
 */
double compute_analytical_gdab(double (*params)[50], const double *arr_z, double *out, size_t n)
{
    const double A = (*params)[0];
    const double H = (*params)[1];
    const double ETA = (*params)[2];

    // Parameter checks
    if (A <= 0.0) {
        fprintf(stderr, "Error: A(%f) <= 0\n", A);
        return -1;
    };

    if (H <= -0.5) {
        fprintf(stderr, "Error: H(%f) <= -1/2\n", H);
        return -1;
    };

    // Precompute scalars
    double V = pow(2.0 * A, 3.0) * M_PI * sqrt(M_PI) * sf_poch(H, 1.5);

    double denom = tgamma(1.5 + H) * 2.0 * M_PI * (A * A);

    double common_prefactor = ETA * ETA;

    // Order of the Bessel function, and G(z -> 0), which is the z -> 0 limit
    // of Gz below (since (u/2)^nu * K_nu(u) -> Gamma(nu)/2).
    const double nu = 0.5 + H;
    const double G0 = (V * V) / (2.0 * M_PI * (A * A) * (1.0 + 2.0 * H));

    // Loop over z
    for (size_t i = 0; i < n; i++) {
        double zi = arr_z[i];

        // Invalid input: negative z
        if (zi < 0.0) {
            fprintf(stderr, "Error: negative value z(%f) found in z_arr\n", zi);
            return -1;
        }

        // z == 0 -> return 0 (Python scalar behavior extended to array)
        if (zi == 0.0) {
            out[i] = 0.0;
            continue;
        }

        double u = zi / A;

        // cyl_bessel_k takes a real order, so integer and non-integer nu use
        // the same call. Do not switch to the Kummer form
        // sqrt(pi) * (2u)^nu * exp(-u) * U(H+1, 2H+2, 2u):
        // boost_hypergeometric_u differences two terms of size ~exp(2u) and so
        // returns nonsense for u >~ 15.
        double KH = bessel_Knu(nu, u);

        // KH alone overflows (as NaN) for small u and large nu, but the
        // product stays bounded and tends to Gamma(nu)/2. Clamping KH to 0
        // instead would give -eta^2 * G0 where the true value is ~0.
        double bessel_term = isfinite(KH) ? pow(u / 2.0, nu) * KH : tgamma(nu) / 2.0;

        // Compute Gz ----------------------------------------------------------
        double Gz = bessel_term * V * V / denom;

        // Final result
        out[i] = common_prefactor * (Gz - G0);
    }

    return 1;
}