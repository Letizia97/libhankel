
#ifndef ANALYTICAL_H
#define ANALYTICAL_H

#include <stddef.h>

double compute_analytical_spheres(
    double (*params)[50], 
    double *arr_z, 
    double *G,
    size_t n
); 

#endif // ANALYTICAL_H