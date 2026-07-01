
#ifndef FORM_FACTORS_H
#define FORM_FACTORS_H

#include "libhankel.h"
#include <stddef.h>


typedef struct {
    double *params;
} form_factor_ctx;


/** 
 * @brief Computes the g_dab form factor. 
 *
 * @param q        value at which to compute the transform
 * @param f_ctx    pointer to struct containing inputs for form factor
 * 
 * @note Params must contain :
 *      -   XI  correlation length (e.g., 10.0)
 *      -   H   Hurst exponent (e.g., 0.5)
 *      -   ETA scattering length density contrast (e.g., 1e-4)
 */
double form_factor_g_dab(double q, void *f_ctx); 



/** 
 * @brief Computes the sphere form factor. 
 *
 * @param q        value at which to compute the transform
 * @param f_ctx    pointer to struct containing inputs for form factor
 * 
 * @note Params must contain :
 *      -   R   radius (e.g., 10.0)
 *      -   ETA scattering contrast (e.g., 1.0)
 */    
double form_factor_sphere(double q, void *f_ctx);



/** 
 * @brief Computes the broad_peak form factor. 
 *
 * @param q        value at which to compute the transform
 * @param f_ctx    pointer to struct containing inputs for form factor
 * 
 * @note Params must contain :
 *      -   I0 forward scattering (e.g., 10e5)
 *      -   XI correlation length (e.g., 1000)
 *      -   Q0 peak position which is related to the d-spacing as q0 = 2pi/d (e.g., 0.01)
 *      -   M  I(Q)=I0/(1+(|q-q0|*xi)^m)^p (e.g., 2)
 *      -   P  I(Q)=I0/(1+(|q-q0|*xi)^m)^p (e.g., 2)
 */  
double form_factor_broad_peak(double q, void *f_ctx);

form_factor_f get_form_factor_by_name(const char *name);

#endif // FORM_FACTORS_H
