
#ifndef QWE_KEY
#define QWE_KEY
#include "libhankel.h"

// sasfit strategy 12
double qwe_Key(int nu, form_factor_f f, double x, void *f_params, double *output, int n_max_iters,
               double rtol, double atol);

#endif // QWE_KEY