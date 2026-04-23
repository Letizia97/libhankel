
#ifndef INTDE_H
#define INTDE_H


// Needed for sasfit strategy 0 (DE QUADRATURE)
void sasfit_intdeini(int lenaw, double tiny, double eps, double *aw);
void sasfit_intde(double (*f)(double, void *), double a, double b, double *aw, 
    double *i, double *err, void *f_params);
void sasfit_intdeoini(int lenaw, double tiny, double eps, double *aw);
void sasfit_intdeo(double (*f)(double, void *), double a, double omega, double *aw, 
    double *i, double *err, void *f_params);

#endif // INTDE_H