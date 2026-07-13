// #include "libhankel.h"

// #include <math.h>
// #include <stdio.h>
// #include <gsl/gsl_math.h>
// #include <gsl/gsl_sf.h>
// #include <gsl/gsl_sf_bessel.h>

// /*
// This file contains functions corresponding to strategies 2-4 in SASfit.
// Specifically:
// - HANKEL_FBT0 2
// - HANKEL_FBT1 3
// - HANKEL_FBT2 4
// */

// // Struct containing both the function to hankel-transform
// // and the function parameters
// struct f_context {
//     double (*f)(double, double (*)[50]);
//     double (*f_params)[50];
// };

// // This wrapping of the function f is needed to be
// // able to multiply by x before passing f to the transform
// double wrapped_f(double x, void *ctx) {
//     struct f_context *c = ctx;
//     return x * c->f(x, c->f_params);
// }

// /**
//  * @brief Computes Hankel transform using FBT.
//  * @note Corresponds to strategies 2,3,4 in SASfit or HANKEL_FBT0, etc
//  * @note Requires knowledge of the form factor function, as one of the
//  parameters
//  *       passed (f_max) is a starting guess for maximum in form factor.
//  *       Performance is highly dependent on params.
//  *
//  * @param nu         order of bessel function - must be 0 or 1
//  * @param f          pointer to form factor function
//  * @param x          value at which to compute the transform
//  * @param f_params   params for form factor
//  * @param n_eval     integer indicating number of function evaluations
//  (N_ogata in SASfit)
//  * @param f_max      float indicating starting guess for maximum in form
//  factor (h_ogata in SASfit)
//  * @param n_method   FBT method to use among 0,1,2 (modified, unmodified,
//  fixed h Ogata)
//  */
// double hankel_transform_FBT(
//     int nu,
//     double (*f)(double, double (*)[50]),
//     double x,
//     double (*f_params)[50],
//     double *output,
//     int n_method,
//     int n_eval,
//     double f_max) {

//     struct f_context ctx = { f, f_params };
//     *output = compute_hankel_FBT(nu, wrapped_f, &ctx, x, n_method, n_eval,
//     f_max); return 0;
// }
