
#ifndef TANH_SINH_QUAD_H
#define TANH_SINH_QUAD_H

// used in sasfit_integrate_ctm
double TanhSinhQuad(double (*f)(double, void *), void *par, double a, double b, int n, double eps,
                    double *err);

#endif // TANH_SINH_QUAD_H