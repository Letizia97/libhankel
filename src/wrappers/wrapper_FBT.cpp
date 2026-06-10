// bessel.cpp
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <chrono>
#include "FBT.h"
extern "C"
{
    #include "include/libhankel.h"
}


double compute_hankel_FBT(
    int nu, 
    form_factor_f intern_fct,
    void *ctx,
    double x, 
    int n_method, 
    int n_eval, 
    double f_max
){
    /*
    Computes Hankel transform using function from external FBT library.

    Receives:
        - nu:          order of Bessel function
        - intern_fct:  function to use for computing the transform, meant to receive x
                       and *ctx (see below for what these are)
        - ctx:         struct containing both the function f to hankel-transform and 
                       the function parameters
        - x:           point at which to compute the transform
        - n_method:    integer among 0,1,2 for modified, unmodified, fixed h Ogata
        - n_eval:      constant corresponding to number of function calls 
        - f_max:       constant corresponding to initial guess where function has maximum   
    */

    // This is a bit of a hack, but it's needed to avoid issues when calling ogata functions
    gsl_set_error_handler_off();

    FBT FBT_instance = FBT(nu, n_method, n_eval, f_max);
	return FBT_instance.fbt(std::bind(intern_fct, std::placeholders::_1, ctx), x);
} 