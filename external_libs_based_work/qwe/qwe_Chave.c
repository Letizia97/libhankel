
/*
Copyright (c) 2018, Alan Chave
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the distribution
    * Neither the name of the Woods Hole Oceanographic Institution nor the names
      of its contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/


#include "external_libs_based_work/qwe/qwe.h"

#include <float.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <stdbool.h>

#include "src/utils/sasfit_integrate.h"


double bessel_j_zero(int nzero, double order) {
    /*
    This function is called Zeroj in SASfit.
    Computes an approximation to the nth positive zero of the Bessel function
    of the first kind, for the special cases where the order is 0 or 1
      If nzero < 7 the algorithm in J. Comp. Phys., 42, 403-405, 1981 is used
      If nzero >= 7, McMahon's asymptotic expansion is used

    Receives:
        nzero  index of the zero, integer >= 1
        order  order of the bessel function, must be 0 or 1    
    */
    double j_zero, beta, t, b;

    if (nzero < 1) {
        fprintf(stderr, 
            "nzero, that is the index of the zero to be computed,"
            " must be greater than or equal to 1\n"
        );
        return -1;
    }

    if (nzero < 7) {
        beta = (nzero + order/2.0 - 0.25) * M_PI;

        if (order == 0) {
            t = 0.0682894897349453
                + 0.131420807470708 * beta * beta
                + 0.0245988241803681 * gsl_pow_4(beta)
                + 0.000813005721543268 * gsl_pow_6(beta);
            b = beta
                + 1.16837242570470 * gsl_pow_3(beta)
                + 0.200991122197811 * gsl_pow_5(beta) 
                + 0.00650404577261471 * gsl_pow_7(beta);
            j_zero = beta + t/b;
            return j_zero;

        } else if (order == 1) {
            t = - 0.362804405737084 
                + 0.120341279038597 * gsl_pow_3(beta)
                + 0.0439454547101171 * gsl_pow_4(beta) 
                + 0.00159340088474713 * gsl_pow_6(beta);
            b = beta 
                - 0.325641790801361 * gsl_pow_3(beta) 
                - 0.117453445968927 * gsl_pow_5(beta) 
                - 0.424906902601794 * gsl_pow_7(beta);
            j_zero = beta + t/b;
            return j_zero;
        }

    } else {

        if (order == 0) {
            beta = (nzero - 0.25) * M_PI;
            j_zero = beta 
                   + 1.0 / (8.0 * beta) 
                   - 124.0 / (1536.0 * gsl_pow_3(beta)) 
                   + 120928.0 / (491520.0 * gsl_pow_5(beta)) 
                   - 401743168.0 / (220200960.0 * gsl_pow_7(beta));
            return j_zero;

        } else if (order == 1) {
            beta = (nzero + 0.25) * M_PI;
            j_zero = beta 
                   - 3.0 / (8.0 * beta) 
                   + 36.0 / (1536.0 * gsl_pow_3(beta)) 
                   - 113184.0 / (491520.0 * gsl_pow_5(beta)) 
                   + 1951209.0 / (220200960.0 * gsl_pow_7(beta));
            return j_zero;
        }
    }
    // FIXME: I think this should be removed completely, and 
    // just support order 0 or 1 (and so add an if statement
    // at the start for this)
    return gsl_sf_bessel_zero_Jnu(order, nzero);
}


double pade_sum(double *s, int n) {
    /*
    Computes sum from 1 to n of s(i) using Padé approximant implemented with
    continued fraction expansion; see Z. Naturforschung, 33a, 402-417, 1978.

    Receives:
        s    series of values to be summed, may be complex   //FIXME : says complex, but only doubles allowed
        n    end of summation
    Outputs:
        sum_cf   sum of the series    
    */

    double *D;     // intermediate “modified moments” (linear combinations of the input terms)
    double *d;     // continued‑fraction coefficients
    double *x, *t; // temporary workspaces used to update polynomial/CF terms
    double sum_cf;
    int i, k, L;

    if (n < 1) {
        fprintf(stderr, "n passed as input to pade_sum must be >= 1.\n");
        return -1;
    }

    D = calloc(n + 1, sizeof(double));
    x = calloc(n + 1, sizeof(double));
    d = calloc(n + 1, sizeof(double));
    t = calloc(n + 1, sizeof(double));


    if (!D || !d || !x || !t) {
        // Allocation failed, free any successful allocations
        free(D); free(d); free(x); free(t);
        fprintf(stderr, "Failed to allocate internal variables in function pade_sum.\n");
        return -1;    
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
        return -1;  
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

        // FIXME: is this fine or should there be a fallback?
        if (D[i-1] == 0) {
            fprintf(stderr, "Division by 0 encountered in pade_sum!\n");
            return -1;  
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




double qwe_Chave(
    double nu, 
    double (*f)(double, double (*)[50]), 
    double r, 
    void *fparams, 
    int n_iters, 
    double rtol, 
    double atol
) {
    /*
    Computes Hankel transform integral from 0 to inf J_sub_order(x*r)*f(x) by integration
    between zero crossings of the Bessel function followed by summation
    using Pade approximant to speed up convergence

    Input variables
        nu       order of the Bessel function, either 0 or 1
        f        function to compute kernel called as f(x,fparams)
        r        x where to compute the Hankel transform
        fparams  input params for f
        n_iters  max number of partial integral intervals
        rtol     relative error
        atol     absolute error
 
    Output variable
        res     computed integral
    */
    double res, a, b, last_res, *s, req_accuracy;
    int nzero;
	bool converged;
	hankel_inputs inputs;
	converged = false;
    last_res = 0;

    // FIXME: add a guard to check nu is 0 or 1 ?

	inputs.function = f;
	inputs.other_inputs[0] = nu;
	inputs.other_inputs[1] = r;
    inputs.fparams = fparams;

    b = bessel_j_zero(1, nu) / r*rtol;    
    s = calloc(n_iters+1, sizeof(double));

    if (!s) {
        // Allocation failed
        fprintf(stderr, "Failed to allocate internal variables in function sasfit_HankelChave.\n");
        return -1;    
    }

    // FIXME: are we sure this is safe?
    for (nzero=1; nzero<=n_iters; nzero++) {     //upper limit is arbitrary and should never be reached
        a = b;
        b = bessel_j_zero(nzero, nu) / r;
        s[nzero] = sasfit_integrate_ctm(a, b, &FrJnu, &inputs, 10000, atol, rtol);
        res = pade_sum(s, nzero);

        req_accuracy = rtol*fabs(res) + atol;
        if (fabs(res-last_res) <= req_accuracy) {
            converged = true;
            return res;
        }
        last_res = res;
    }
    free(s);

    if (!converged) {
        fprintf("HankelChave algorithm did not converge after maximum allowed intervals: %d\n",n_iters);
        return -1;
    };

    return res;
}
