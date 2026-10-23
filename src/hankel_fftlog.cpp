#include "../external_libs/fft_log/fft_log.h"
#include <math.h>

extern "C" {
#include "libhankel.h"
#include "utils/validate_x.h"
}

extern "C" int hankel_transform_fftlog(int nu, form_factor_f f, const double *x, size_t len_x,
                             void *f_ctx, double *output) {
    int status = validate_x_array(x, len_x);
    if (status != 0) return status;

    for (size_t j = 0; j < len_x; j++) {
        const double v = sasfit_fftlog((double)nu, f, x[j], f_ctx);
        if (isnan(v)) return -1;
        output[j] = v;
    }
    return 0;
}
