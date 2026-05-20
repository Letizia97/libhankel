#include "src/utils/sf_functions.h"
#include <math.h>
#include "libhankel.h"

double sf_poch(double a, double x) {
    return tgamma(a + x) / tgamma(a);
}