#include "libhankel.h"
#include "../external_libs/qwe_Chave/qwe_Chave.h"
#include "../external_libs/qwe_Key/qwe_Key.h"

// Standard library headers
#include <float.h>
#include <math.h>
#include <stdio.h>


/*
This file contains functions corresponding to strategies 12 and 13 in SASfit. 
They have been grouped together as these are QWE transforms, that is Quadrature With Extrapolation.
Specifically:
- HANKEL_QWE 12
- HANKEL_CHAVE 13
*/



double hankel_transform_QWE_Key(
    int nu, 
    form_factor_f f, 
    const double x,
    void *f_ctx,
    double * output,
    int n_eval, 
    double eps_rel
) {

    int status;
    status = qwe_Key(
        nu, 
        f, 
        x, 
        f_ctx,
        output,
        lround(n_eval),
        eps_rel*10,
        DBL_MIN);
        
    return status;
}


double hankel_transform_QWE_Chave(
    int nu, 
    form_factor_f f, 
    const double x,
    void *f_ctx,
    double * output,
    int n_eval, 
    double eps_rel
) {
    int status;
    status = qwe_Chave(
        nu, 
        f,
        x, 
        f_ctx,
        output,
        lround(n_eval),
        eps_rel*10,
        DBL_MIN);

    return status;
}
