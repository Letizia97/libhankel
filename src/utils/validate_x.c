#include "src/utils/validate_x.h"

// Standard library headers
#include <math.h>
#include <stddef.h>
#include <stdio.h>

int validate_x(double x) {
    if (!isfinite(x) || x <= 0) {
        fprintf(stderr, "Error: x must be finite and greater than zero\n");
        return -12;
    }
    return 0;
}

int validate_x_array(const double *x, size_t len_x) {
    for (size_t j = 0; j < len_x; j++) {
        int status = validate_x(x[j]);
        if (status != 0) {
            return status;
        }
    }
    return 0;
}
