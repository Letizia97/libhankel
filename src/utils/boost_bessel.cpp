// bessel_wrapper.cpp
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/hypergeometric_1F1.hpp>
#include <cmath>
#include <exception>

// These wrappers are declared extern "C" and called from C code. Boost.Math
// uses its default (throwing) policy, so invalid arguments (e.g. non-positive
// x, k < 1, gamma poles) raise C++ exceptions. Letting an exception propagate
// through C stack frames is undefined behavior and in practice calls
// std::terminate(). We therefore catch here and return NaN, which callers can
// detect with isnan() and which propagates harmlessly through the numerics.

extern "C" {

double bessel_Jnu(double nu, double x) {
    try {
        return boost::math::cyl_bessel_j(nu, x);
    } catch (const std::exception &) {
        return NAN;
    }
}

double bessel_Jnu_zero(double nu, int k) {
    try {
        return boost::math::cyl_bessel_j_zero(nu, k);
    } catch (const std::exception &) {
        return NAN;
    }
}

double bessel_Knu(double nu, double x) {
    try {
        return boost::math::cyl_bessel_k(nu, x);
    } catch (const std::exception &) {
        return NAN;
    }
}

double boost_hypergeometric_u(double a, double b, double x) {
    using boost::math::hypergeometric_1F1;

    try {
        return std::tgamma(1.0 - b) / std::tgamma(a - b + 1.0) * hypergeometric_1F1(a, b, x) +
               std::tgamma(b - 1.0) / std::tgamma(a) * std::pow(x, 1.0 - b) *
                   hypergeometric_1F1(a - b + 1.0, 2.0 - b, x);
    } catch (const std::exception &) {
        return NAN;
    }
}
}
