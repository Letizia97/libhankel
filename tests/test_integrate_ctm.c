#include <libhankel.h>
#include "unity.h"
#include "utils/unity_config.h"

// Standard library headers
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

// Project / local headers
#include "form_factors.h"
#include "src/utils/sasfit_integrate.h"
#include "utils/test_utils.h"


#define MAX_COLS 25
#define MAX_ROWS 6

double nu;
double z;

form_factor_ctx ctx_spheres;
form_factor_ctx ctx_gdab;

hankel_inputs inputs;
double G0_spheres;
double G0_gdab;

// For running the tests
typedef struct {
    double actual_spheres;
    double expected_spheres;
    double actual_gdab;
    double expected_gdab;
} TestContext;

TestContext ctx = {
    .expected_spheres = 394781.721568501,
    .actual_spheres = 0, 
    .expected_gdab = 0.0315827340834857,
    .actual_gdab = 0, 
};

// Needed for computing G0
double spheres_ff_at_0(double q, hankel_inputs *params){
    return q * form_factor_sphere(q, params);
}

double gdab_ff_at_0(double q, hankel_inputs *params){
    return q * form_factor_g_dab(q, params);
}


void setUp(void) {
    /*
    Function that runs before tests, for a shared set-up.
    Anything computed inside it is available in tests
    because previously declared outside the function.        
    */
    
    // set params
    double params_spheres[] = {10.0, 1.0};
    ctx_spheres.params = params_spheres;

    double params_gdab[] = {10.0, 0.5, 1e-4};
    ctx_gdab.params = params_gdab;
    
    nu = 0;

    // set inputs
    inputs.function = form_factor_sphere;
    inputs.f_params = params_spheres;
	inputs.other_inputs[0] = nu;
	inputs.other_inputs[1] = 0;

    G0_spheres = sasfit_integrate_ctm(
        0,
        INFINITY,
        spheres_ff_at_0,
        &inputs,
        10000,
        0.001,
        1e-20
    );
    ctx.actual_spheres = G0_spheres;

    inputs.function = form_factor_g_dab;
    inputs.f_params=params_gdab;
    G0_gdab = sasfit_integrate_ctm(
        0,
        INFINITY,
        gdab_ff_at_0,
        &inputs,
        10000,
        0.001,
        1e-20
    );
    ctx.actual_gdab = G0_gdab;

}


void tearDown(void) {}

void test_integrate_ctm_spheres(void) {
    /*
    Regression test for sasfit_integrate_ctm on spheres
    */
    TEST_ASSERT_EQUAL(ctx.actual_spheres, ctx.expected_spheres);
}

void test_integrate_ctm_gdab(void) {
    /*
    Regression test for sasfit_integrate_ctm on gdab
    */
    TEST_ASSERT_EQUAL(ctx.actual_gdab, ctx.expected_gdab);
}


int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_integrate_ctm_spheres);
    RUN_TEST(test_integrate_ctm_gdab);
    return UNITY_END();
}
