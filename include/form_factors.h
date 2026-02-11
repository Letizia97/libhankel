
#ifndef FORM_FACTORS_H
#define FORM_FACTORS_H


/** 
 * @brief Computes the g_dab form factor. 
 *
 * @param q          value at which to compute the Hankel t
 * @param params     pointer to an array of params for the function
 * 
 * @note Params must contain :
 *      -   XI
 *      -   H
 *      -   ETA
 */
double form_factor_g_dab(double q, double (*params)[50]); 



/** 
 * @brief Computes the sphere form factor. 
 *
 * @param q          value at which to compute the Hankel t
 * @param params     pointer to an array of params for the function
 * 
 * @note Params must contain :
 *      -   R
 *      -   ETA
 */    
double form_factor_sphere(double q, double (*params)[50]);



/** 
 * @brief Computes the broad_peak form factor. 
 *
 * @param q          value at which to compute the Hankel t
 * @param params     pointer to an array of params for the function
 * 
 * @note Params must contain :
 *      -   I0
 *      -   XI
 *      -   Q0
 *      -   M
 *      -   P
 */ 
double form_factor_broad_peak(double q, double (*params)[50]);


#endif // FORM_FACTORS_H
