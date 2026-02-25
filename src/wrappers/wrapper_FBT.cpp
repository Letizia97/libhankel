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


double compute_hankel_FBT(int nu, double (*intern_fct)(double, double (*)[50]), double x, double (*f_params)[50], int n_method, int N, double h){
    /*
    Computes Hankel transform using function from external FBT library.

    Receives:
        - x:    point at which to compute the transform
        - (*intern_fct)(double, void *):    function to use for computing the transform
        - * f_params:    parameters for the function
        - nu: order of Bessel function
        - n_method: integer among 0,1,2 for : modified, unmodified, fixed h Ogata
        - N: constant corresponding to number of function calls 
        - h: constant corresponding to initial guess where function has maximum   
    */

    // This is a bit of a hack, but it's needed to avoid issues when calling ogata functions
    gsl_set_error_handler_off();

    FBT FBT_instance = FBT(nu, n_method, N, h);
	return FBT_instance.fbt(std::bind(intern_fct, std::placeholders::_1, f_params), x);
} 