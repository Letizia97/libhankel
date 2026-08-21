
#ifndef QWE_KEY
#define QWE_KEY
#include "libhankel.h"

#include <stddef.h>

// sasfit strategy 12
int qwe_Key(int nu, form_factor_f f, const double *x, size_t len_x, void *f_params, double *output,
            int n_max_iters, double rtol, double atol);

#endif // QWE_KEY