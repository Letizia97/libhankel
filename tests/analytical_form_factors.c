#include "tests/analytical_form_factors.h"

#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include "include/libhankel.h"



double compute_analytical_spheres(double (*params)[50], double *arr_z, double *G, size_t n) {
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
            //printf("Entered block A\n");

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
    } else {
        return 0;      // added
    }
}
