#include "test_common.h"  

#include <stdio.h>
#include <math.h>

int arrays_close(double *actual,
                        double *expected,
                        size_t n,
                        double tol)
{
    for (size_t i = 0; i < n; ++i) {
        if (fabs(actual[i] - expected[i]) > tol) {
            return 0;
        }
    }
    return 1;
}
