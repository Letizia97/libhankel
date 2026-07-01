#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <stdbool.h>
#include "utils/unity_config.h"
#include "unity.h"
#include <string.h>

#include "libhankel.h"
#include "form_factors.h"
#include "src/utils/sasfit_integrate.h"
#include "utils/test_utils.h" 

#define ARRAY_LEN 19

double nu;
double z;

form_factor_ctx ctx_spheres;
form_factor_ctx ctx_gdab;
form_factor_ctx ctx_broad_peak;

double r_array_spheres[ARRAY_LEN];
double r_array_gdab[ARRAY_LEN];
double r_array_broad_peak[ARRAY_LEN];

hankel_inputs inputs;
double Gr[ARRAY_LEN];



void setUp(void) {
    
    nu = 0;

    // set params
    double params_spheres[] = {10.0, 1.0};
    ctx_spheres.params = params_spheres;

    double params_gdab[] = {10.0, 0.5, 1e-4};
    ctx_gdab.params = params_gdab;

    double params_broad_peak[] = {10e5, 1000, 0.01, 2, 2};
    ctx_broad_peak.params = params_broad_peak;

    strategy_params strategy_params_general = {
        .n_eval = 250, 
        .eps_rel = 1e-9
    };

    strategy_params strategy_params_b_peak = {
        .n_eval = 150, 
        .eps_rel = 1e-9
    };

    // setup the x (or r) array
    double r_array_spheres[ARRAY_LEN] = {
        1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0,  10.0,
        11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0
    };

    double r_array_gdab[ARRAY_LEN] = {
        15.,          18.54166667,  22.08333333,  25.625,       29.16666667,
        32.70833333,  36.25,        39.79166667,  43.33333333,  46.875,
        50.41666667,  53.95833333,  57.5,         61.04166667,  64.58333333,
        68.125,       71.66666667,  75.20833333,  78.75,             
    };

    double r_array_broad_peak[ARRAY_LEN] = {
        30.,   49.,   69.,   88.,   108.,  127.,  147.,  167.,  186.,  206.,  225., 
        225.,  245.,  265.,  284.,  304.,  323.,  343.,  362.,       
    };

    // COMPUTATIONS 
    hankel_transform(
        nu, 
        form_factor_sphere, 
        r_array_spheres, 
        ARRAY_LEN,
        (void *)&ctx_spheres,
        ctx.actual_spheres, 
        "QWE_Key", 
        strategy_params_general
    );

    hankel_transform(
        nu, 
        form_factor_g_dab, 
        r_array_gdab, 
        ARRAY_LEN,
        (void *)&ctx_gdab,
        ctx.actual_gdab, 
        "QWE_Key", 
        strategy_params_general
    );

    hankel_transform(
        nu, 
        form_factor_broad_peak, 
        r_array_broad_peak, 
        ARRAY_LEN,
        (void *)&ctx_broad_peak,
        ctx.actual_broad_peak, 
        "QWE_Key", 
        strategy_params_b_peak
    );
}
