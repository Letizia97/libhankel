
#ifndef QWE_H
#define QWE_H

// used in sasfit_integrate_ctm
double TanhSinhQuad(
	double (*f)(double, void *), 
	void *par, 
	double a, 
	double b, 
	int n, 
	double eps, 
	double *err
);

// sasfit strategy 12
double qwe_Key(
	double nu, 
	double (*f)(double, double (*)[50]), 
	double x, 
	void *fparams, 
	int nIntervalsMax, 
	double rtol, 
	double atol
); 

// sasfit strategy 13
double qwe_Chave(
    double nu, 
    double (*f)(double, double (*)[50]), 
    double r, 
    void *fparams, 
    int n_iters, 
    double rtol, 
    double atol
);


#endif // QWE_H