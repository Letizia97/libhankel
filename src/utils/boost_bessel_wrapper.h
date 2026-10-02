// bessel_wrapper.h
#ifndef BESSEL_WRAPPER_H
#define BESSEL_WRAPPER_H

// Tests call these functions directly, so they must be exported from the DLL on Windows.
#include "libhankel_export.h"

LIBHANKEL_API double bessel_Jnu(double nu, double x);
LIBHANKEL_API double bessel_Jnu_zero(double nu, int k);
LIBHANKEL_API double bessel_Knu(double nu, double x);
LIBHANKEL_API double boost_hypergeometric_u(double a, double b, double x);
#endif
