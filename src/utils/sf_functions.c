#include "src/utils/sf_functions.h"
#include "libhankel.h"
#include <math.h>

double sf_poch(double a, double x) { return tgamma(a + x) / tgamma(a); }