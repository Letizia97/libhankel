
#include "../external_libs/DE-quadrature/intde.h"
#include "libhankel.h"

// Standard library headers
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Project / local headers
#include "../src/utils/boost_bessel_wrapper.h"

/*
This file contains only 2 Hankel strategies, both DE quadrature algorithms,
corresponding to:
    - HANKEL_OOURA_DEO , i.e. strategy 0 in SASfit
    - OGATA_2005 , i.e. strategy 1 in SASfit
*/

/** Length of the node/weight table handed to the Ooura integrators.
 *  Ooura requires lenaw > 1000; 8000 is his suggested value for IEEE double. */
#define DE_WORKSPACE_LEN 4000

/**
 * @brief Everything the integrand needs, threaded through the Ooura
 *        integrators as an opaque `void *`.
 */
typedef struct {
    form_factor_f f; /**< function to integrate */
    void *f_ctx;     /**< parameters for the supplied function */
    int nu;          /**< order of the Bessel function */
    double Q;        /**< radial Fourier variable, i.e. conj wavenumber to radius.
                          This is the `x` argument of the public entry points. */
} hankel_integrand_ctx;

/**
 * @brief Auxiliary function that computes the value of the
 *        Hankel‑transform integrand at the current radius r.
 *        Computes a product between the radius r, the Bessel
 *        function of the first kind of order nu, and a function
 *        supplied through ctx.
 *
 * @param r    radius
 * @param ctx  pointer to a @ref hankel_integrand_ctx containing nu, Q and the
 *             function to integrate
 */
static double hankel_integrand(double r, void *ctx) {
    hankel_integrand_ctx *integrand_ctx = (hankel_integrand_ctx *)ctx;
    if (r == 0)
        return 0;

    int nu = integrand_ctx->nu;
    double Q = integrand_ctx->Q;

    double bessel = jn(nu, Q * r);
    double fval = integrand_ctx->f(r, integrand_ctx->f_ctx);

    return r * bessel * fval;
}

/**
 * @brief Auxiliary function computing the
 *        double‑exponential (DE) / tanh–sinh transform.
 */
static double DEtransform(double t) {
    double s = sinh(t);
    double a = M_PI_2 * s;
    return t * tanh(a);
}

/**
 * @brief Auxiliary function that computes the derivative of a
 *        “double‑exponential (DE) change of variables”, that is
 *        commonly used in numerical integration (quadrature).
 *        Specifically, this is the tanh–sinh (double exponential)
 *        transformation.
 */
static double deriv_DEtransform(double t) {
    double sh = sinh(t);
    double ch = cosh(t);
    double A = M_PI_2 * sh;
    double secH = 1.0 / cosh(A); // sech(A)

    return M_PI_2 * t * ch * (secH * secH) + tanh(A);
}

int hankel_transform_DE_Ooura(int nu, form_factor_f f, const double x, void *f_ctx, double *output,
                              int n_eval, double eps_rel) {

    double res0, err0, res, err;

    /* Both the finite and the oscillatory piece are set up in units of 1/x
     * (see scaled_zero below), so a non-positive or non-finite x would divide
     * by zero and hand NaN to the integrators instead of failing. */
    if (!(x > 0)) {
        fprintf(stderr, "Error: x must be finite and greater than zero\n");
        return -12;
    }

    hankel_integrand_ctx integrand_ctx;
    integrand_ctx.f_ctx = f_ctx;
    integrand_ctx.f = f;
    integrand_ctx.nu = nu;
    integrand_ctx.Q = x;

    // chooses which zero index to request, caps it at 10
    int zero_index = n_eval < 10 ? n_eval : 10;

    // compute zero through bessel function and scale it
    double scaled_zero = bessel_Jnu_zero(nu, zero_index) / x;

    // allocates an array of doubles and returns a pointer to it
    double *workspace = malloc(DE_WORKSPACE_LEN * sizeof *workspace);
    if (workspace == NULL) {
        fprintf(stderr, "Failed to allocate internal variables "
                        "in function hankel_transform_DE_Ooura.\n");
        return -3;
    }

    // precompute nodes & weights for DE integration on a finite interval [a, b]
    sasfit_intdeini(DE_WORKSPACE_LEN, DBL_MIN, eps_rel, workspace);
    // compute integral using DE quadrature with weights created by
    // sasfit_intdeini
    sasfit_intde(&hankel_integrand, 0, scaled_zero, workspace, &res0, &err0, &integrand_ctx);
    // precompute nodes/weights for oscillatory integrals
    // e.g. f(x) cos(omega x) , over [a, inf]
    sasfit_intdeoini(DE_WORKSPACE_LEN, DBL_MIN, eps_rel, workspace);
    // evaluate oscillatory integrals using the table built by intdeoini.
    sasfit_intdeo(&hankel_integrand, scaled_zero, integrand_ctx.Q, workspace, &res, &err,
                  &integrand_ctx);

    free(workspace);
    res += res0;

    // estimate of the numerical integration error for the computed integral
    // FIXME: need to figure out what to do with this
    err += err0;

    *output = res;
    return 0;
}

/* Note on the argument names, which follow libhankel.h and SASfit rather than
 * Ogata's paper: `x` is the radial Fourier variable (Q below), and `f_max`
 * (`h_ogata` in SASfit) is a starting guess for the maximum of
 * q * formFactor(q).  It enters the quadrature below in the position of
 * Ogata's step size h, so the scale of the integrand's peak is what sets the
 * node spacing. */
int hankel_transform_DE_Ogata(int nu, form_factor_f f, const double x, void *f_ctx, double *output,
                              int n_eval, double f_max) {

    double sum;
    sum = 0.0;

    /* The quadrature nodes are y_k / x and the result carries a 1 / x^2, so a
     * non-positive or non-finite x would silently produce Inf or NaN. */
    if (!(x > 0)) {
        fprintf(stderr, "Error: x must be finite and greater than zero\n");
        return -12;
    }

    for (int i = 1; i <= n_eval; i++) {

        /* ---- Get Bessel zero α_{ν,i} scaled by π ---- */
        double zero_i = bessel_Jnu_zero(nu, i);
        double zero_scaled = zero_i / M_PI;

        /* ---- Apply DE transform & its derivative ---- */
        double t = f_max * zero_scaled;
        double phi = DEtransform(t);
        double phi_prime = deriv_DEtransform(t); /* Jacobian */

        /* ---- Map DE node into actual integration node y_k ---- */
        double y_k = phi * (M_PI / f_max);

        /* ---- Precompute Bessel factors ---- */
        double Jnu_yk = jn(nu, y_k);
        double Jnu1_zero = jn(nu + 1, zero_i);

        /* ---- Quadrature weight for α_{ν,i} ---- */
        double denom = M_PI * Jnu1_zero;
        double weight = 2.0 / ((denom * denom) * zero_scaled);

        /* ---- Evaluate integrand at scaled location ---- */
        double f_val = (*f)(y_k / x, f_ctx);

        /* ---- Assemble quadrature contribution ---- */
        double term = weight * y_k * f_val * Jnu_yk * phi_prime;

        sum += term;
    }

    // Apply normalization from the change of variables t = x * r
    double scaled_sum = (M_PI / (x * x)) * sum;

    // Parity correction for integer Bessel order (H_-ν = (-1)^ν H_ν).
    // Only bites for negative odd nu; hankel_transform() restricts nu to
    // {0, 1}, so this is a no-op unless this function is called directly.
    double parity = (nu < 0 && nu % 2 != 0) ? -1.0 : 1.0;

    // Final result
    double res = scaled_sum * parity;

    *output = res;
    return 0;
}
