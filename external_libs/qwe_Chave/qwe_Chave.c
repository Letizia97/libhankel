#include "external_libs/qwe_Chave/qwe_Chave.h"
#include "libhankel.h"

// Standard library headers
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

// Project / local headers
#include "src/utils/pow_functions.h"
#include "src/utils/sasfit_integrate.h"


/** 
 * @brief This function is called Zeroj in SASfit.
 * @note Computes an approximation to the nth positive zero of the Bessel function
 *       of the first kind, for the special cases where the order is 0 or 1
 *       If nzero < 7 the algorithm in J. Comp. Phys., 42, 403-405, 1981 is used
 *       If nzero >= 7, McMahon's asymptotic expansion is used
 * 
 * @param nzero  index of the zero, integer >= 1
 * @param order  order of the bessel function, must be 0 or 1    
 */
double bessel_j_zero(int nzero, double order) {

    double j_zero, beta, t, b;

    if (!(order==0 || order==1)) {
        fprintf(stderr, "order needs to be 0 or 1 in order to use bessel_j_zero\n");
        return -1;
    }

    if (nzero < 1) {
        fprintf(stderr, 
            "nzero, that is the index of the zero to be computed, "
            "must be greater than or equal to 1\n"
        );
        return -6;
    }

    if (nzero < 7) {
        beta = (nzero + order/2.0 - 0.25) * M_PI;

        if (order == 0) {
            t = 0.0682894897349453
                + 0.131420807470708 * beta * beta
                + 0.0245988241803681 * pow4(beta)
                + 0.000813005721543268 * pow6(beta);
            b = beta
                + 1.16837242570470 * pow3(beta)
                + 0.200991122197811 * pow5(beta) 
                + 0.00650404577261471 * pow7(beta);
            j_zero = beta + t/b;
            return j_zero;

        } else if (order == 1) {
            t = - 0.362804405737084 
                + 0.120341279038597 * pow3(beta)
                + 0.0439454547101171 * pow4(beta) 
                + 0.00159340088474713 * pow6(beta);
            b = beta 
                - 0.325641790801361 * pow3(beta) 
                - 0.117453445968927 * pow5(beta) 
                - 0.424906902601794 * pow7(beta);
            j_zero = beta + t/b;
            return j_zero;
        }

    } else {

        if (order == 0) {
            beta = (nzero - 0.25) * M_PI;
            j_zero = beta 
                   + 1.0 / (8.0 * beta) 
                   - 124.0 / (1536.0 * pow3(beta)) 
                   + 120928.0 / (491520.0 * pow5(beta)) 
                   - 401743168.0 / (220200960.0 * pow7(beta));
            return j_zero;

        } else if (order == 1) {
            beta = (nzero + 0.25) * M_PI;
            j_zero = beta 
                   - 3.0 / (8.0 * beta) 
                   + 36.0 / (1536.0 * pow3(beta)) 
                   - 113184.0 / (491520.0 * pow5(beta)) 
                   + 1951209.0 / (220200960.0 * pow7(beta));
            return j_zero;
        }
    }
    return 0;
}

/** 
 * @brief Computes sum from 1 to n of s(i) using Padé approximant implemented with
 *        continued fraction expansion; see Z. Naturforschung, 33a, 402-417, 1978.
 * 
 * @param s    series of values to be summed
 * @param n    end of summation
 */
double pade_sum(double *s, int n) {
    double *D;     // intermediate “modified moments” (linear combinations of the inputs)
    double *d;     // continued‑fraction coefficients
    double *x, *t; // temporary workspaces used to update polynomial/CF terms
    double sum_cf;
    int i, k, L;

    if (n < 1) {
        fprintf(stderr, "n passed as input to pade_sum must be >= 1.\n");
        return -7;
    }

    D = calloc(n + 1, sizeof(double));
    x = calloc(n + 1, sizeof(double));
    d = calloc(n + 1, sizeof(double));
    t = calloc(n + 1, sizeof(double));


    if (!D || !d || !x || !t) {
        // Allocation failed, free any successful allocations
        free(D); free(d); free(x); free(t);
        fprintf(stderr, 
            "Failed to allocate internal variables "
            "in function pade_sum.\n"
        );
        return -3;    
    }

    D[1] = s[1];
    d[1] = D[1];
    if (n == 1) {
        sum_cf = d[1];
        goto cleanup_and_exit;
    }
    D[2] = s[2];

    if (D[1] == 0) {
        fprintf(stderr, 
            "Division by 0 encountered in pade_sum!"
            " s[1] must be different from 0. \n"
        );
        return -5;  
    }

    d[2] = -D[2]/D[1];
    if (n == 2) {
        sum_cf = d[1] / (1 + d[2]);
        goto cleanup_and_exit;
    }
    for (i=3; i<=n; i++) {
        L = 2 * lround(floor((i-1.)/2.));
        //update x vector
        for (k=L; k>=4; k=k-2) {
            x[k] = x[k-1] + d[i-1] * x[k-2];
        }
        x[2] = x[1] + d[i - 1];
        // interchange odd and even parts
        for (k=1; k<=L-1; k=k+2) {
            t[k] = x[k];
            x[k] = x[k+1];
            x[k+1] = t[k];
        }

        // compute cf coefficient
        D[i] = s[i];
        for (k=0; k<=L/2-1; k++) {
            D[i] = D[i] + s[i-1-k] * x[1 + 2*k];
        }
        //    D[i] = s[i] + s[i-1:-1:i-L/2]*x[1:2:L-1];

        if (D[i-1] == 0) {
            fprintf(stderr, "Division by 0 encountered in pade_sum!\n");
            return -5;  
        }
        d[i] = -D[i]/D[i-1];   
    }
    //evaluate continued fraction
    sum_cf = 1;
    for (k=n; k>=2; k--) {
        sum_cf = 1 + d[k]/sum_cf;
    }

    sum_cf = d[1]/sum_cf;
cleanup_and_exit:
    free(D);
    free(d);
    free(x);
    free(t);
    return sum_cf;
}


/** 
 * @brief Computes Hankel transform integral from 0 to inf J_sub_order(x*r)*f(x)
 *        by integration between zero crossings of the Bessel function followed 
 *        by summation using Pade approximant to speed up convergence
 * 
 * @param nu           order of the Bessel function, either 0 or 1
 * @param f            function to compute kernel called as f(x,f_params)
 * @param r            x where to compute the Hankel transform
 * @param f_params      input params for f
 * @param output       pointer to var containing output from transform 
 * @param n_max_iters  max number of partial integral intervals
 * @param rtol         relative error
 * @param atol         absolute error
 */
double qwe_Chave(
    double nu, 
    form_factor_f f, 
    double r, 
    void *f_params, 
    double *output, 
    int n_max_iters, 
    double rtol, 
    double atol
) {

    double res, a, b, last_res, *s, req_accuracy;
    int nzero;
	bool converged = false;
    last_res = 0;

    if (!(nu==0 || nu==1)) {
        fprintf(stderr, 
            "nu needs to be 0 or 1 in order to use "
            "QWE_Chave\n"
        );
        return -1;
    }

    hankel_inputs inputs;
	inputs.function = f;
	inputs.other_inputs[0] = nu;
	inputs.other_inputs[1] = r;
    inputs.f_params = f_params;

    b = bessel_j_zero(1, nu) / r*rtol;    
    s = calloc(n_max_iters+1, sizeof(double));

    if (!s) {
        // Allocation failed
        fprintf(stderr, 
            "Failed to allocate internal variables in "
            "function sasfit_HankelChave.\n"
        );
        return -3;    
    }

    //upper limit (n_max_iters) is arbitrary and should never be reached
    for (nzero=1; nzero<=n_max_iters; nzero++) {     
        a = b;
        b = bessel_j_zero(nzero, nu) / r;
        s[nzero] = sasfit_integrate_ctm(a, b, &FrJnu, &inputs, 10000, atol, rtol);
        res = pade_sum(s, nzero);

        req_accuracy = rtol*fabs(res) + atol;
        if (fabs(res-last_res) <= req_accuracy) {
            converged = true;
            *output = res;
            return 0;
        }
        last_res = res;
    }
    free(s);

    if (!converged) {
        fprintf(stderr, 
            "QWE_Chave algorithm did not converge "
            "after maximum allowed intervals (%d)\n", n_max_iters
        );
        return -4;
    };

    *output = res;
	return 0;
}
