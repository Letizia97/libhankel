#ifndef UTILS_H
#define UTILS_H

#include <stddef.h>

int arrays_close(
    double *actual,
    double *expected,
    size_t n,
    double tol
);

void start_capture_stderr(void);
void stop_capture_stderr(
    char *buffer, 
    size_t bufsize
);

double* read_values_by_rows(const char *filename,
                            const char *column_name,
                            const double *rows,
                            size_t nrows,
                            double *out);

#endif // UTILS_H