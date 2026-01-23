
#ifndef TEST_UTILS_H
#define TEST_UTILS_H

#include <stddef.h>

static void start_capture_stderr(void);
static void stop_capture_stderr(char *buffer, size_t bufsize);

#endif // TEST_UTILS_H