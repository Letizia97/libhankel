
#ifndef QWE_CHAVE
#define QWE_CHAVE

// sasfit strategy 13
double qwe_Chave(
    double nu, 
    double (*f)(double, double (*)[50]), 
    double r, 
    void *fparams, 
    double *output,
    int n_max_iters, 
    double rtol, 
    double atol
);


#endif // QWE_CHAVE