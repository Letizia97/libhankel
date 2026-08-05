
#include "../external_libs/DE-quadrature/intde.h"
#include "libhankel.h"

// Standard library headers
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

// Project / local headers
#include "../src/utils/boost_bessel_wrapper.h"
#include "../src/utils/validate_x.h"

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

int hankel_transform_DE_Ooura(int nu, form_factor_f f, const double *x, size_t len_x, void *f_ctx,
                              double *output, int n_eval, double eps_rel) {

    /* A non-positive or non-finite x would hand NaN to the integrators. */
    int status = validate_x_array(x, len_x);
    if (status != 0) {
        return status;
    }

    /* Two tables: the finite and oscillatory ones are built from the same
     * inputs but are not interchangeable. Each depends only on eps_rel and is
     * read-only to the integrators, so both are built once for all points. */
    double *aw_finite = malloc(DE_WORKSPACE_LEN * sizeof *aw_finite);
    double *aw_oscillatory = malloc(DE_WORKSPACE_LEN * sizeof *aw_oscillatory);
    if (aw_finite == NULL || aw_oscillatory == NULL) {
        free(aw_finite);
        free(aw_oscillatory);
        fprintf(stderr, "Failed to allocate internal variables "
                        "in function hankel_transform_DE_Ooura.\n");
        return -3;
    }

    // precompute nodes & weights for DE integration on a finite interval [a, b]
    sasfit_intdeini(DE_WORKSPACE_LEN, DBL_MIN, eps_rel, aw_finite);
    // precompute nodes/weights for oscillatory integrals
    // e.g. f(x) cos(omega x) , over [a, inf]
    sasfit_intdeoini(DE_WORKSPACE_LEN, DBL_MIN, eps_rel, aw_oscillatory);

    // chooses which zero index to request, caps it at 10
    const int zero_index = n_eval < 10 ? n_eval : 10;
    // the zero itself does not depend on x; only its scaling below does
    const double zero = bessel_Jnu_zero(nu, zero_index);

    for (size_t j = 0; j < len_x; j++) {
        hankel_integrand_ctx integrand_ctx;
        integrand_ctx.f_ctx = f_ctx;
        integrand_ctx.f = f;
        integrand_ctx.nu = nu;
        integrand_ctx.Q = x[j];

        const double scaled_zero = zero / x[j];

        double res0, err0, res, err;

        // compute integral using DE quadrature with weights created by
        // sasfit_intdeini
        sasfit_intde(&hankel_integrand, 0, scaled_zero, aw_finite, &res0, &err0, &integrand_ctx);
        // evaluate oscillatory integrals using the table built by intdeoini.
        sasfit_intdeo(&hankel_integrand, scaled_zero, integrand_ctx.Q, aw_oscillatory, &res, &err,
                      &integrand_ctx);

        res += res0;

        // estimate of the numerical integration error for the computed integral
        // FIXME: need to figure out what to do with this
        err += err0;

        output[j] = res;
    }

    free(aw_finite);
    free(aw_oscillatory);
    return 0;
}

/* Argument names follow libhankel.h and SASfit rather than Ogata's paper: `x`
 * is the radial Fourier variable (Q below), and `f_max` (`h_ogata` in SASfit)
 * is a starting guess for the maximum of q * formFactor(q). It enters the
 * quadrature in the position of Ogata's step size h, so the scale of the
 * integrand's peak sets the node spacing. */
/**
 * @brief One quadrature node, holding everything that does not depend on the
 *        transform variable x.
 */
typedef struct {
    double y_k;        /**< node; the form factor is sampled at y_k / x */
    double weight_y_k; /**< quadrature weight, premultiplied by y_k */
    double Jnu_yk;     /**< J_nu(y_k) */
    double phi_prime;  /**< Jacobian of the DE transform at this node */
} ogata_node;

/**
 * @brief Builds the @p n_nodes x-independent quadrature nodes.
 *
 * The Bessel zeros are found by a root search and dominate the runtime, which
 * is why this is hoisted out of the loop over x.
 */
static void build_ogata_nodes(ogata_node *nodes, size_t n_nodes, int nu, double f_max) {
    for (size_t k = 0; k < n_nodes; k++) {

        /* ---- Get Bessel zero α_{ν,i} scaled by π ---- */
        double zero_i = bessel_Jnu_zero(nu, (int)k + 1);
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

        nodes[k].y_k = y_k;
        nodes[k].weight_y_k = weight * y_k;
        nodes[k].Jnu_yk = Jnu_yk;
        nodes[k].phi_prime = phi_prime;
    }
}

int hankel_transform_DE_Ogata(int nu, form_factor_f f, const double *x, size_t len_x, void *f_ctx,
                              double *output, int n_eval, double f_max) {

    /* The nodes are y_k / x and the result carries a 1 / x^2, so a
     * non-positive or non-finite x would silently produce Inf or NaN. */
    int status = validate_x_array(x, len_x);
    if (status != 0) {
        return status;
    }

    /* A non-positive n_eval means no nodes and a sum of zero, as the original
     * loop bound gave; this also keeps a negative count out of the malloc. */
    const size_t n_nodes = n_eval > 0 ? (size_t)n_eval : 0;

    ogata_node *nodes = NULL;
    if (n_nodes > 0) {
        nodes = malloc(n_nodes * sizeof *nodes);
        if (nodes == NULL) {
            fprintf(stderr, "Failed to allocate internal variables "
                            "in function hankel_transform_DE_Ogata.\n");
            return -3;
        }
        build_ogata_nodes(nodes, n_nodes, nu, f_max);
    }

    // Parity correction for integer Bessel order (H_-ν = (-1)^ν H_ν).
    // Only bites for negative odd nu; hankel_transform() restricts nu to
    // {0, 1}, so this is a no-op unless this function is called directly.
    const double parity = (nu < 0 && nu % 2 != 0) ? -1.0 : 1.0;

    for (size_t j = 0; j < len_x; j++) {
        double sum = 0.0;

        for (size_t k = 0; k < n_nodes; k++) {
            /* ---- Evaluate integrand at scaled location ---- */
            double f_val = (*f)(nodes[k].y_k / x[j], f_ctx);

            /* ---- Assemble quadrature contribution ---- */
            double term = nodes[k].weight_y_k * f_val * nodes[k].Jnu_yk * nodes[k].phi_prime;

            sum += term;
        }

        // Apply normalization from the change of variables t = x * r
        double scaled_sum = (M_PI / (x[j] * x[j])) * sum;

        output[j] = scaled_sum * parity;
    }

    free(nodes);
    return 0;
}
