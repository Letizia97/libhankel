
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

// #include "include/libhankel.h"

double Zeroj(int nzero,double order) {
    //Computes the approximate zero of bessel function of the first kind of
    //order 0 or 1
    //   If nzero < 7 the algorithm in J. Comp. Phys., 42, 403-405, 1981 is used
    //   If nzero >= 7, McMahon's asymptotic expansion is used

    //Input variables
    //nzero--index of the zero, integer >= 1
    //order--order of the bessel function, must be 0 or 1

    double Zero,beta, t,b;
    if (nzero < 7) {
        beta = (nzero + order/2.0 - 0.25)*M_PI;
        if (order == 0) {
            t = 0.0682894897349453 + 0.131420807470708*beta*beta +
                0.0245988241803681 *gsl_pow_4(beta) +
                0.000813005721543268 *gsl_pow_6(beta);
            b = beta + 1.16837242570470*gsl_pow_3(beta) +
                0.200991122197811*gsl_pow_5(beta) + 0.00650404577261471*gsl_pow_7(beta);
            Zero = beta + t/b;
            return Zero;
        } else if (order == 1) {
            t = -0.362804405737084 + 0.120341279038597*gsl_pow_3(beta)+
                0.0439454547101171*gsl_pow_4(beta) + 0.00159340088474713*gsl_pow_6(beta);
            b = beta - 0.325641790801361*gsl_pow_3(beta) -
                0.117453445968927*gsl_pow_5(beta) - 0.424906902601794*gsl_pow_7(beta);
            Zero = beta + t/b;
            return Zero;
        }
    } else {
        if (order == 0) {
            beta = (nzero - 0.25)*M_PI;
            Zero = beta + 1.0/(8.0*beta) - 124.0/(1536.0*gsl_pow_3(beta)) +
                120928.0/(491520.0*gsl_pow_5(beta)) -
                401743168.0/(220200960.0*gsl_pow_7(beta));
            return Zero;
        } else if (order == 1) {
            beta = (nzero + 0.25)*M_PI;
            Zero = beta - 3.0/(8.0*beta) + 36.0/(1536.0*gsl_pow_3(beta)) -
                113184.0/(491520.0*gsl_pow_5(beta)) +
                1951209.0/(220200960.0*gsl_pow_7(beta));
            return Zero;
        }
    }
    return gsl_sf_bessel_zero_Jnu(order,nzero);
}


double Padesum(double *s, int n) {
    //Computes sum from 1 to n of s(i) using Pade approximant implemented with
    //continued fraction expansion; see Z. Naturforschung, 33a, 402-417, 1978.

    //Input variable
    //s--series of values to be summed, may be complex
    //Output variable
    //Cf--sum of the series

    double *D, *d, *x, *t;
    double Cf;
    int i,k,L;
    D = calloc(n+1,sizeof(double));
    x = calloc(n+1,sizeof(double));
    d = calloc(n+1,sizeof(double));
    t = calloc(n+1,sizeof(double));
    D[1] = s[1];
    d[1] = D[1];
    if (n == 1) {
        Cf = d[1];
        goto exitPadesum;
    }
    D[2] = s[2];
    d[2] = -D[2]/D[1];
    if (n == 2) {
        Cf = d[1]/(1 + d[2]);
        goto exitPadesum;
    }
    for (i=3;i<=n;i++) {
        L = 2*lround(floor((i-1.)/2.));
        //update x vector
        for (k=L;k>=4;k=k-2) {
            x[k] = x[k-1] + d[i-1]*x[k-2];
        }
        x[2] = x[1] + d[i - 1];
        // interchange odd and even parts
        for (k=1;k<=L-1;k=k+2) {
            t[k] = x[k];
            x[k] = x[k+1];
            x[k+1] = t[k];
        }

        // compute cf coefficient
        D[i] = s[i];
        for (k=0;k<=L/2-1;k++) {
            D[i] = D[i] + s[i-1-k]*x[1+2*k];
        }
        //    D[i] = s[i] + s[i-1:-1:i-L/2]*x[1:2:L-1];
        d[i] = -D[i]/D[i-1];
    }
    //evaluate continued fraction
    Cf = 1;
    for (k=n;k>=2;k--) {
        Cf = 1 + d[k]/Cf;
    }

    Cf = d[1]/Cf;
exitPadesum:
    free(D);
    free(d);
    free(x);
    free(t);
    return Cf;
}




double sasfit_HankelChave(double order, double (*f)(double, double (*)[50]), double r, void *fparams, int nIntervalsMax, double rerr, double aerr) {
    /*
    Computes Hankel transform integral from 0 to inf J_sub_order(x*r)*f(x) by integration
    between zero crossings of the Bessel function followed by summation
    using Pade approximant to speed up convergence

    Input variables
        order   order of the Bessel function, either 0 or 1
        r       argument of the Hankel transform
        f       function to compute kernel called as f(x,fparams)
        rerr    relative error
        aerr    absolute error

    Output variable
        Sum     computed integral
    */
    double Sum,a,b,last,*s;
    int nzero;
	bool converged;
	hankel_inputs inputs;
	converged = false;
	inputs.function = f;
	inputs.other_inputs[0] = order;
	inputs.other_inputs[1] = r;
    inputs.fparams=fparams;
    b=0;
    b=Zeroj(1,order)/r*rerr;
    last=0;
    s=calloc(nIntervalsMax+1,sizeof(double));
    for (nzero = 1;nzero<=nIntervalsMax;nzero++) {     //upper limit is arbitrary and should never be reached
        a = b;
        b = Zeroj(nzero,order)/r;
        s[nzero] = sasfit_integrate_ctm(a, b,&FrJnu,&inputs, 10000, aerr, rerr);
        Sum = Padesum(s,nzero);
        if (fabs(Sum - last) <= rerr*fabs(Sum) + aerr) {
            converged = true;
            return Sum;
        }
        last=Sum;
    }
    free(s);

    // FIXME: proper error handling
    if (!converged) printf("HankelChave algorithm did not converge after maximum allowed intervals: %d\n",nIntervalsMax);
    return Sum;
}
