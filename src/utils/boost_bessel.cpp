// bessel_wrapper.cpp
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/hypergeometric_1F1.hpp>

extern "C" {

double bessel_Jnu(double nu, double x) { return boost::math::cyl_bessel_j(nu, x); }

double bessel_Jnu_zero(double nu, int k) { return boost::math::cyl_bessel_j_zero(nu, k); }

double bessel_Knu(double nu, double x) { return boost::math::cyl_bessel_k(nu, x); }

double boost_hypergeometric_u(double a, double b, double x) {
    using boost::math::hypergeometric_1F1;

    return std::tgamma(1.0 - b) / std::tgamma(a - b + 1.0) * hypergeometric_1F1(a, b, x) +
           std::tgamma(b - 1.0) / std::tgamma(a) * std::pow(x, 1.0 - b) *
               hypergeometric_1F1(a - b + 1.0, 2.0 - b, x);
}
}
