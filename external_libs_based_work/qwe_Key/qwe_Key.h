
#ifndef QWE_KEY
#define QWE_KEY

// sasfit strategy 12
double qwe_Key(
	double nu, 
	double (*f)(double, double (*)[50]), 
	double x, 
	void *fparams, 
    double *output,
	int n_max_iters, 
	double rtol, 
	double atol
); 


#endif // QWE_KEY