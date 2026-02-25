/*
%------------------------------------------------------------------%
% Copyright (c) 2012 by the Society of Exploration Geophysicists.  %
% For more information, go to http://software.seg.org/2012/0003 .  %
% You must read and accept usage terms at:                         %
% http://software.seg.org/disclaimer.txt before use.               %
%------------------------------------------------------------------%

function T = qwe(rtol,atol,nIntervalsMax,func,varargin)
%
% QWE - Quadrature with extrapolation
%
% Computes an infinite integral using the partial sum of quadrature terms
% accelerated by sequence extrapolation using the Shanks transformation
% implemented with Wynn's epsilon algorithm.
%
% Usage:
%
% T = qwe(rtol,atol,nIntervalsMax,func,varargin)
%
% other_inputs:
%
% rtol,atol     - Relative and absolute tolerances for the extrapolation.
%                 Extrapolation for integral T stops when:
%                 absolute_error < rtol*abs(T) + atol
% nIntervalsMax - Maximum number of partial integral intervals
% func,varargin - Handle to the function that returns the
%                 quadarature sum for this interval. varargin
%                 are passed through to func to assist in the kernel eval.
%                 func has the form: func(interval,varargin), where the
%                 interger interval is an index to the precomputed absiccae
%                 and weights stored in the calling routine. func can
%                 return multiple related kernels, each of which will be
%                 integrated and extrapolated.
%
% Outputs:
%
% T - structure with fields described below. If func returns multiple
% related kernels (nKernels), then the T is an array of structures, ie
% T(1:nKernels).  The kernel integrations and subsequent transforms are
% applied to each kernel in parallel, but are stopped uniquely for each
% when tolerance has  been met.  For example, func might return the J0 and
% J1 kernels of electromagnetic fields Ey (nKernels=2), or from multiple
% components (nKernels=2*3 for Ex,Ey,Ez). This routine would then return
% T(1).n, T(2).n etc with the fields shown below:
%
% T.n             - Final order of the transform to reach rtol or atol.
% T.extrap(1:T.n) - Result from extrapolation at each order.  The final
%                   best answer is: T.extrap(T.n)
% T.absErr(1:T.n) - Estimated absolute error for each order n.
% T.relErr(1:T.n) - Estimated relative error for each order n.
%
%
% Written by:
%
% Kerry Key
% Scripps Institution of Oceanography
%
% History:
%
% June 8, 2011 - Moved quadrature outside.
% May 2011     - Implemented.
%
%--------------------------------------------------------------------------

%
% Adjustable parameters:
%
    nDelay = 1;  % must be >= 1, but 1 is usually sufficient.

%
% First compute nDelay partial integrals before starting
% the Shanks transformation iterations:
%
    prev = 0;
    for i = 1:nDelay
        prev    = prev + func(i,varargin{:}); % Call to func returns quadrature sum
    end

%
% Then initialize the T structure for the extrapolation results:
%
    nKernels = length(prev);  % func from nDelay loop above returns nKernels results, so now we
                              % know how many related kernels are being handled

    nTerms   = nIntervalsMax-nDelay-1;  % maximum number of terms in extrapolation

    rmin     = realmin;

    for i = 1:nKernels  % Preallocating arrays at maximum size for speed:
        T(i).S      = zeros(nTerms,1); % working array used for the recursion coefficients for the Epsilon algorithm
        T(i).extrap = zeros(nTerms,1); % extrapolated result for each order of the expansion
        T(i).relErr = zeros(nTerms,1); % relative error for each order
        T(i).absErr = zeros(nTerms,1); % absolute error for each order
        converged   = false(nKernels,1);
    end

%
% The extrapolation transformation loop:
%
    for i = nDelay+1:nTerms

        %
        % Step 1: Compute Guass quadrature of this interval:
        %
        f = func(i,varargin{:});

        %
        % Step 2: Compute Shanks transformation for each kernel function:
        %
        for j = 1:nKernels

            if  converged(j)  % skip this kernel since its done
                continue;
            end

            %
            % Insert components into the T structure:
            %
            n           = i-nDelay;
            T(j).n      = n;               % order of the expansion
            T(j).S(n+1) = T(j).S(n)+f(j);  % working array for transformation

            %
            % Compute the Shanks transform using the Epsilon algorithm:
            % Structured after Weniger (1989, p26)
            %
            aux2 = 0;
            for k = n+1:-1:2
                aux1 = aux2;
                aux2 = T(j).S(k-1);
                ddff = T(j).S(k) - aux2;
                if abs(ddff) < rmin
                    T(j).S(k-1) = realmax;
                else
                    T(j).S(k-1) = aux1 + 1/ddff;
                end
            end

            %
            % The extrapolated result plus the prev integration term:
            %
            T(j).extrap(n) = T(j).S(mod(n,2)+1)+prev(j);

            %
            % Step 3: Analyze for convergence:
            %
            if n > 1
                T(j).absErr(n) = abs( T(j).extrap(n) - T(j).extrap(n-1));
                T(j).relErr(n) = T(j).absErr(n) / abs(T(j).extrap(n)) ;
                converged(j)   = T(j).relErr(n) < rtol + atol/abs(T(j).extrap(n));
            end

        end % loop over nKernels

        if all(converged)
            break;
        end

    end % extrapolation loop

    % Clean up the T structure arrays:
    for j = 1:nKernels
        n = T(j).n;
        T(j).extrap = T(j).extrap(1:n);
        T(j).relErr = T(j).relErr(1:n);
        T(j).absErr = T(j).absErr(1:n);
    end

    % Goodbye my friend, returning to the calling function

end % function qwe
*/

#include "external_libs_based_work/qwe_Key/qwe_Key.h"

#include <float.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <stdbool.h>

#include "src/utils/sasfit_integrate.h"

/** 
 * @brief Computes Hankel transform integral using strategy 13 from SASfit
 * 
 * @param nu           order of the Bessel function, either 0 or 1
 * @param f            function to compute kernel called as f(x,fparams)
 * @param r            x where to compute the Hankel transform
 * @param fparams      input params for f
 * @param output       pointer to var containing output from transform 
 * @param n_max_iters  max number of partial integral intervals
 * @param rtol         relative error
 * @param atol         absolute error
 */
double qwe_Key(
    double nu, 
    double (*f)(double, double (*)[50]), 
    double x, 
    void *fparams, 
    double *output,
    int n_max_iters,
    double rtol, 
    double atol
) {

	int idx_of_zero = 1; // index of the zero to compute, must be >= 1, but 1 is usually sufficient
	int n_terms;  // maximum number of terms in extrapolation
	double last_res = 0;  //latest res in iters
	double *S; // array used for the recursion coefficients for the Epsilon algorithm
    double *extrap;  // extrapolated result for each order of the expansion
    double *rel_err, *abs_err; // rel and abs errors 
    double r_min, r_max; // start and end of integration range 
    int i, k, n;
    double f_i, aux2, aux1, diff, res;
	bool converged = false;

    if (!(nu==0 || nu==1)) {
        fprintf(stderr, 
            "nu needs to be 0 or 1 in order to use "
            "the selected strategy\n"
        );
        return -1;
    }

    hankel_inputs inputs;
	inputs.function = f;
	inputs.other_inputs[0] = nu;
	inputs.other_inputs[1] = x;
    inputs.fparams=fparams;

	r_max = gsl_sf_bessel_zero_Jnu(nu, idx_of_zero) / inputs.other_inputs[1];
	r_min = r_max * (rtol / 10);

    // First compute idx_of_zero partial integrals before starting the Shanks transformation iters
	last_res = last_res + sasfit_integrate_ctm(r_min, r_max, &FrJnu, &inputs, 10000, atol, rtol);
	n_terms = n_max_iters - idx_of_zero - 1; 

	S = calloc(n_terms+1,sizeof(double));
	extrap = calloc(n_terms+1,sizeof(double));
	rel_err = calloc(n_terms+1,sizeof(double));
	abs_err = calloc(n_terms+1,sizeof(double));

    if (!S || !extrap || !rel_err || !abs_err) {
        // Allocation failed, free any successful allocations
        free(S); free(extrap); free(rel_err); free(abs_err);
        fprintf(stderr, 
            "Failed to allocate internal variables "
            "in function pade_sum.\n"
        );
        return -1;    
    }

    for (i=idx_of_zero+1; i<=n_terms; i++) {
        r_min = r_max;
        r_max = gsl_sf_bessel_zero_Jnu(nu, i) / inputs.other_inputs[1];

        // compute Guass quadrature of this interval
        f_i = sasfit_integrate_ctm(r_min, r_max, &FrJnu, &inputs, 10000, atol, rtol);

        n = i - idx_of_zero; // order of the expansion
        S[n+1] = S[n] + f_i;

        // Compute the Shanks transform using the Epsilon algorithm:
        // Structured after Weniger (1989, p26)
        aux2 = 0.0;
        for (k=n+1; k>=2; k--) {
            aux1 = aux2;
            aux2 = S[k-1];
            diff = S[k] - aux2;
            if (fabs(diff) < DBL_MIN) {
                S[k-1] = DBL_MAX;
            } else {
                S[k-1] = aux1 + 1./diff;
            }
        }

        extrap[n] = S[(n % 2)+1] + last_res;
        res = extrap[n];

        if (n > 1) {
            abs_err[n] = fabs( extrap[n] - extrap[n-1]);
            rel_err[n] = abs_err[n] / fabs(extrap[n]) ;
            if (rel_err[n] < rtol + atol/fabs(extrap[n]))
                converged = true;
        }
        if (converged) break;
    }
	free(S);
	free(extrap);
	free(rel_err);
	free(abs_err);

    if (!converged) {
        printf(
            "qwe_Key algorithm did not converge "
            "after maximum allowed intervals: %d\n",n_max_iters
        );
        return -1;
    };

    *output = res;
	return 0;
}

