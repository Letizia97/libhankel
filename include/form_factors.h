
#ifndef FORM_FACTORS_H
#define FORM_FACTORS_H

#include "libhankel.h"
#include <stddef.h>

/**
 * @struct form_factor_ctx
 * @brief Context structure for form factor callbacks.
 *
 * This structure is used to pass user-defined data to functions of type
 * @ref form_factor_f. The contents of the context are managed by the user
 * and interpreted by the callback implementation.
 */
typedef struct {

    /**
     * @brief Pointer to user-defined parameter array.
     *
     * This array stores the data required by the form factor
     * callback. For indication of what this should contain,
     * please refer to the documentation for each form factor.
     */
    double *params;
} form_factor_ctx;

/**
 * @brief Computes the g_dab form factor.
 *
 * @param q        value at which to compute the transform
 * @param f_ctx    pointer to struct containing inputs for form factor
 *                 ``f_ctx`` must contain a pointer to an array with the
 *                 following parameters, in this order:
 *                  -   \f$ \xi \f$: correlation length (e.g., 10.0)
 *                  -   H : Hurst exponent (e.g., 0.5)
 *                  -   ETA : scattering length density contrast (e.g., 1e-4)
 */
double form_factor_g_dab(double q, void *f_ctx);

/**
 * @brief Computes the sphere form factor.
 *
 * @param q        value at which to compute the transform
 * @param f_ctx    pointer to struct containing inputs for form factor
 *                 ``f_ctx`` must contain a pointer to an array with the
 *                 following parameters, in this order:
 *                  -   R : radius (e.g., 10.0)
 *                  -   ETA : scattering contrast (e.g., 1.0)
 */
double form_factor_sphere(double q, void *f_ctx);

/**
 * @brief Computes the broad_peak form factor.
 *
 * @param q        value at which to compute the transform
 * @param f_ctx    pointer to struct containing inputs for form factor.
 *                 ``f_ctx`` must contain a pointer to an array of parameters.
 *                 The parameters, which appear in the formula
 *                 \f[ I(Q) = \frac{I_0}{\left(1 + (|q - q_0| \xi)^m \right)^p}
 * \f] need to be in the following order:
 *                   - \f$ I_0 \f$ : Forward scattering intensity (e.g. 1e5)
 *                   - \f$ \xi \f$ : Correlation length (e.g. 1000)
 *                   - \f$ q_0 \f$: Peak position, related to d-spacing:
 *                     \f$ q_0 = \frac{2\pi}{d} \f$ (e.g. 0.01)
 *                   - \f$ m \f$ : Exponent \f$m\f$ (e.g. 2)
 *                   - \f$ p \f$ : Exponent \f$p\f$ (e.g. 2)
 */
double form_factor_broad_peak(double q, void *f_ctx);

form_factor_f get_form_factor_by_name(const char *name);

#endif // FORM_FACTORS_H
