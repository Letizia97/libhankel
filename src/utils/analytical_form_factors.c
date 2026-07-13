#include "src/utils/analytical_form_factors.h"

#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>

#include "libhankel.h"

#include "../src/utils/sf_functions.h"

double compute_analytical_spheres(double (*params)[50], const double *arr_z, double *G, size_t n) {
    /*
    Unnormalized SESANS correlation function for a homogeneous sphere of radius R.
    Parameters
    ----------
    normalized : bool, optional
        If True, return G(z)/G(0). Otherwise return unnormalized G(z).
    log_form : {'A','B'}, optional
        Choose the logarithm argument:
        'A' -> ln( 1 + xi / sqrt(1 - xi**2) )
        'B' -> ln( (1 + sqrt(1 - xi**2)) / xi )
        Default 'B' tends to be numerically safer near δ→2R.

    Notes
    -----
    - The function is defined for 0 <= δ <= 2R. For δ > 2R, the correlation is 0.
    - Numerical safeguards are applied near xi -> 1 to avoid overflow.
    */

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

double compute_analytical_gdab(double (*params)[50], const double *arr_z, double *out, size_t n) {
    const double A = (*params)[0];
    const double H = (*params)[1];
    const double ETA = (*params)[2];

    // Parameter checks
    if (A < 0.0) {
        fprintf(stderr, "Error: A(%f) < 0\n", A);
        return -1;
    };

    if (H <= -0.5) {
        fprintf(stderr, "Error: H(%f) <= -1/2\n", H);
        return -1;
    };

    // Precompute scalars
    double V = pow(2.0 * A, 3.0) * M_PI * sqrt(M_PI) * sf_poch(H, 1.5);

    double denom = gsl_sf_gamma(1.5 + H) * 2.0 * M_PI * (A * A);

    double common_prefactor = ETA * ETA;

    // Loop over z
    for (size_t i = 0; i < n; i++) {
        double zi = arr_z[i];

        // Invalid input: negative z
        if (zi < 0.0) {
            fprintf(stderr, "Negative value found in z_arr");
            return -1;
        }

        // z == 0 -> return 0 (Python scalar behavior extended to array)
        if (zi == 0.0) {
            out[i] = 0.0;
            continue;
        }

        double u = zi / A;

        double G0 = (V * V) / (2.0 * M_PI * (A * A) * (1.0 + 2.0 * H));

        // Determine if (H + 0.5) is an integer -------------------------------
        bool is_integer_case = fabs((0.5 + H) - round(0.5 + H)) < 1e-12;

        double KH;

        if (is_integer_case) {
            int n_int = (int)round(H + 0.5);

            // GSL Bessel K_v: here v = n_int
            KH = gsl_sf_bessel_Knu((double)n_int, u);
        } else {
            // KH = sqrt(pi) * (2u)^(H + 0.5) * exp(-u) * U(H+1, 2H+2, 2u)
            double power = pow(2.0 * u, H + 0.5);
            double U = gsl_sf_hyperg_U(H + 1.0, 2.0 * H + 2.0, 2.0 * u);

            KH = sqrt(M_PI) * power * exp(-u) * U;
        }

        // Handle invalid KH ---------------------------------------------------
        if (!isfinite(KH)) {
            KH = 0.0;
        }

        // Compute Gz ----------------------------------------------------------
        double Gz = KH * V * V * pow(u / 2.0, 0.5 + H) / denom;

        // Final result
        out[i] = common_prefactor * (Gz - G0);
    }

    return 1;
}