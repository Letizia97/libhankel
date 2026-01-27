#include "utils.h"  

#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#include "unity_config.h"
#include "unity.h"

int arrays_close(double *actual,
                        double *expected,
                        size_t n,
                        double tol)
{
    for (size_t i = 0; i < n; ++i) {
        if (fabs(actual[i] - expected[i]) > tol) {
            return 0;
        }
    }
    return 1;
}

static int saved_stderr = -1;
static int pipefd[2] = {-1, -1};

void start_capture_stderr(void)
{
    fflush(stderr);
    saved_stderr = dup(STDERR_FILENO);
    TEST_ASSERT_NOT_EQUAL(-1, saved_stderr);

    TEST_ASSERT_EQUAL(0, pipe(pipefd));
    // Redirect stderr to pipe write end
    TEST_ASSERT_NOT_EQUAL(-1, dup2(pipefd[1], STDERR_FILENO));
    close(pipefd[1]); // not needed by this process anymore
}

void stop_capture_stderr(char *buffer, size_t bufsize)
{
    // Restore stderr
    fflush(stderr);
    dup2(saved_stderr, STDERR_FILENO);
    close(saved_stderr);

    // Read what was captured
    ssize_t n = read(pipefd[0], buffer, bufsize - 1);
    if (n < 0) n = 0;
    buffer[n] = '\0';
    close(pipefd[0]);
}
