#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <stddef.h>

typedef enum { INTERP_LINEAR, INTERP_PCHIP, INTERP_MAKIMA } interp_method_t;

typedef struct interp interp_t;

interp_t *interp_create(const double *x, const double *y, size_t n, interp_method_t method);

double interp_eval(interp_t *interp, double x);

void interp_destroy(interp_t *interp);

#endif