// clang-format off
#include "utils/unity_config.h"
#include "unity.h"
// clang-format on
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "form_factors.h"
#include "libhankel.h"
#include "src/utils/sasfit_integrate.h"
#include "utils/test_utils.h"

#define ARRAY_LEN 25

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

// For running the tests
typedef struct {
    double actual_spheres[ARRAY_LEN];
    double expected_spheres[ARRAY_LEN];
    double actual_gdab[ARRAY_LEN];
    double expected_gdab[ARRAY_LEN];
    double actual_broad_peak[ARRAY_LEN];
    double expected_broad_peak[ARRAY_LEN];
} TestContext;

TestContext ctx = {
    .expected_spheres =
        {
            387507.478876507,   371194.865559785,  349050.878110533,   322866.912359191,
            293957.890774645,   263387.284336425,  232057.309283767,   200753.920990995,
            170170.357980906,   140920.013857726,  113542.663215219,   88506.3597290707,
            66205.9256771091,   46958.3386835007,  30994.6269460026,   18447.4103664975,
            9331.51981974001,   3512.03954686726,  643.54697816898,    -0.419462003739866,
            0.671565816397894,  0.408622932460936, -0.572610815496937, 0.444648150293901,
            -0.730794091442533,

        },
    .actual_spheres = {0},
    .expected_gdab =
        {
            0.0131409968968235,   0.0099409587529311,   0.00744536154206769,  0.00553395596369854,
            0.00408864031395487,  0.00300616660224412,  0.00220142925163886,  0.00160668986388701,
            0.00116925997763928,  0.000848815153527063, 0.000614860095841608, 0.000444542888287165,
            0.00032086242301512,  0.000231243932588246, 0.000166431099098308, 0.00011963763470055,
            8.59047980625098e-05, 6.16207867381733e-05, 4.41601904858169e-05, 3.16199665583075e-05,
            2.26226721933264e-05, 1.61732894051951e-05, 1.15543192494369e-05, 8.24883138794311e-06,
            5.88501976636725e-06,
        },
    .actual_gdab = {0},
    .expected_broad_peak = {15.3478,  14.7545,  13.8463,  12.7398,  11.3539,  9.86578,   8.1652,
                            6.37775,  4.64776,  2.84404,  1.19288,  1.19288,  -0.433731, -1.904,
                            -3.12269, -4.19016, -4.98079, -5.56571, -5.88294, -5.97027,  -5.81729,
                            -5.4678,  -4.91068, -4.23042, -3.39016},
    .actual_broad_peak = {0},
};

void setUp(void) {
    /*
    Function that runs before tests, for a shared set-up.
    Anything computed inside it is available in tests
    because previously declared outside the function.
    */
    nu = 0;

    // set params
    double params_spheres[] = {10.0, 1.0};
    ctx_spheres.params = params_spheres;

    double params_gdab[] = {10.0, 0.5, 1e-4};
    ctx_gdab.params = params_gdab;

    double params_broad_peak[] = {10e5, 1000, 0.01, 2, 2};
    ctx_broad_peak.params = params_broad_peak;

    // setup the x (or r) array
    double r_array_spheres[ARRAY_LEN] = {1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0,
                                         10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0,
                                         19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0};

    double r_array_gdab[ARRAY_LEN] = {
        15.,         18.54166667, 22.08333333, 25.625,      29.16666667, 32.70833333, 36.25,
        39.79166667, 43.33333333, 46.875,      50.41666667, 53.95833333, 57.5,        61.04166667,
        64.58333333, 68.125,      71.66666667, 75.20833333, 78.75,       82.29166667, 85.83333333,
        89.375,      92.91666667, 96.45833333, 100.,
    };

    double r_array_broad_peak[ARRAY_LEN] = {
        30.,  49.,  69.,  88.,  108., 127., 147., 167., 186., 206., 225., 225., 245.,
        265., 284., 304., 323., 343., 362., 382., 402., 421., 441., 460., 480.,
    };

    // COMPUTATIONS
    // printf("de ogata, Gr-G0  \n");
    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        z = r_array_spheres[i];
        hankel_transform_DE_Ooura(nu, form_factor_sphere, z, (void *)&ctx_spheres, &Gr[i], 250,
                                  1e-9);
        // printf("%.15g, ", (Gr[i]));
        ctx.actual_spheres[i] = Gr[i];
    }

    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        z = r_array_gdab[i];
        hankel_transform_DE_Ooura(nu, form_factor_g_dab, z, (void *)&ctx_gdab, &Gr[i], 250, 1e-9);
        ctx.actual_gdab[i] = Gr[i];
    }

    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        z = r_array_broad_peak[i];
        hankel_transform_DE_Ooura(nu, form_factor_broad_peak, z, (void *)&ctx_broad_peak, &Gr[i],
                                  250, 1e-9);
        ctx.actual_broad_peak[i] = Gr[i];
    }
}

void tearDown(void) {}

void test_hankel_DE_Ooura_regression_spheres(void) {
    /*
    Regression test for QWE on spheres.
    */
    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        TEST_ASSERT_DOUBLE_WITHIN(1e-4, ctx.expected_spheres[i], ctx.actual_spheres[i]);
    }
}

void test_hankel_DE_Ooura_regression_gdab(void) {
    /*
    Regression test for QWE on gdab.
    */
    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        TEST_ASSERT_DOUBLE_WITHIN(1e-4, ctx.expected_gdab[i], ctx.actual_gdab[i]);
    }
}

void test_hankel_DE_Ooura_regression_broad_peak(void) {
    /*
    Regression test for QWE on broad peak.
    */
    for (size_t i = 0; i < ARRAY_LEN; ++i) {
        TEST_ASSERT_DOUBLE_WITHIN(1e-4, ctx.expected_broad_peak[i], ctx.actual_broad_peak[i]);
    }
}

int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_hankel_DE_Ooura_regression_spheres);
    RUN_TEST(test_hankel_DE_Ooura_regression_gdab);
    RUN_TEST(test_hankel_DE_Ooura_regression_broad_peak);
    return UNITY_END();
}
