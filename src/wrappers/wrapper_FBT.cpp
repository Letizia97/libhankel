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


double compute_hankel_FBT(int nu, double (*intern_fct)(double, double (*)[50]), double x, double (*fparams)[50], int n_method, int N, double h){
    /*
    Computes Hankel transform using function from external FBT library.

    Receives:
        - x:    point at which to compute the transform
        - (*intern_fct)(double, void *):    function to use for computing the transform
        - * fparams:    parameters for the function
        - nu: order of Bessel function
        - n_method: integer between 0 and 2
        - N: constant controlling ....
        - h: constant controlling ...
    */
    FBT FBT_instance = FBT(nu, N, h, n_method);
    // return FBT_instance.fbt(std::_Bind_helper<false, double (*&)(double, double (*)[50]), const std::_Placeholder<1>&, void*&>::type, x);
	return FBT_instance.fbt(std::bind(intern_fct, std::placeholders::_1, fparams),x);
} 