// clang-format off
#include "utils/unity_config.h"
#include "unity.h"
// clang-format on
#include <math.h>

/*
Prototypes for the extern "C" Boost.Math wrappers defined in
src/utils/boost_bessel.cpp. These are exercised with invalid arguments to make
sure a C++ exception raised inside Boost is caught at the boundary and turned
into a NaN return value, rather than unwinding through C frames (undefined
behavior / std::terminate).
*/
double bessel_Jnu(double nu, double x);
double bessel_Jnu_zero(double nu, int k);
double bessel_Knu(double nu, double x);

void setUp(void) {}
void tearDown(void) {}

void test_bessel_Knu_returns_nan_for_nonpositive_x(void) {
    /*
    cyl_bessel_k requires x > 0. With Boost's default throwing policy this
    would raise a domain_error; the wrapper must catch it and return NaN.
    */
    TEST_ASSERT_TRUE_MESSAGE(isnan(bessel_Knu(0.0, -1.0)),
                             "bessel_Knu should return NaN for x <= 0");
}

void test_bessel_Jnu_zero_returns_nan_for_invalid_k(void) {
    /*
    The zero index k must be >= 1. An out-of-range k throws in Boost and must
    be caught and reported as NaN.
    */
    TEST_ASSERT_TRUE_MESSAGE(isnan(bessel_Jnu_zero(0.0, 0)),
                             "bessel_Jnu_zero should return NaN for k < 1");
}

void test_bessel_wrappers_valid_for_good_input(void) {
    /*
    Sanity check: valid inputs must still produce finite results (the catch
    block must not swallow legitimate values).
    */
    TEST_ASSERT_FALSE_MESSAGE(isnan(bessel_Knu(0.0, 1.0)), "bessel_Knu(0, 1) should be finite");
    TEST_ASSERT_FALSE_MESSAGE(isnan(bessel_Jnu(0.0, 1.0)), "bessel_Jnu(0, 1) should be finite");
    TEST_ASSERT_FALSE_MESSAGE(isnan(bessel_Jnu_zero(0.0, 1)),
                              "bessel_Jnu_zero(0, 1) should be finite");
}

int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_bessel_Knu_returns_nan_for_nonpositive_x);
    RUN_TEST(test_bessel_Jnu_zero_returns_nan_for_invalid_k);
    RUN_TEST(test_bessel_wrappers_valid_for_good_input);
    return UNITY_END();
}
