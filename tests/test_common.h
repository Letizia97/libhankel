#ifndef TEST_COMMON_H
#define TEST_COMMON_H

#include <stddef.h>

int arrays_close(double *actual,
                        double *expected,
                        size_t n,
                        double tol);

#endif // TEST_COMMON_H