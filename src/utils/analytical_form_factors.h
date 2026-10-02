
#ifndef ANALYTICAL_H
#define ANALYTICAL_H

#include <stddef.h>
// Tests call these functions directly, so they must be exported from the DLL on Windows.
#include "libhankel_export.h"

LIBHANKEL_API double compute_analytical_spheres(double (*params)[50], const double *arr_z, double *G, size_t n);

LIBHANKEL_API double compute_analytical_gdab(double (*params)[50], const double *arr_z, double *out, size_t n);

#endif // ANALYTICAL_H