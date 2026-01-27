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

#endif // UTILS_H