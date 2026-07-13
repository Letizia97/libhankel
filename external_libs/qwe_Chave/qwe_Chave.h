
#ifndef QWE_CHAVE
#define QWE_CHAVE
#include "libhankel.h"

// sasfit strategy 13
double qwe_Chave(int nu, form_factor_f f, double r, void *f_params, double *output, int n_max_iters,
                 double rtol, double atol);

#endif // QWE_CHAVE