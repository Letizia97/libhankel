
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

// Needed for sasfit strategy 12
double sasfit_qwe(
	double nu, 
	double (*f)(double, double (*)[50]), 
	double x, 
	void *fparams, 
	int nIntervalsMax, 
	double rtol, 
	double atol
); 

// Needed for sasfit strategy 13
double sasfit_HankelChave(
	double order, 
	double (*f)(double, double (*)[50]), 
	double r, 
	void *fparams, 
	int nIntervalsMax, 
	double rerr, 
	double aerr
);


#endif // QWE_H