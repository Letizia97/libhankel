// bessel_wrapper.cpp
#include <boost/math/special_functions/bessel.hpp>

extern "C" {

double bessel_Jnu(double nu, double x) {
    return boost::math::cyl_bessel_j(nu, x);
}

double bessel_Jnu_zero(double nu, int k) {
    return boost::math::cyl_bessel_j_zero(nu, k);
}

}