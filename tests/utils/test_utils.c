#include "test_utils.h"  

// External / third-party libraries (Unity test framework)
#include "unity.h"
#include "unity_config.h"

// Standard library headers
#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>



int arrays_close(
    double *actual,
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


/**
 * Reads values from a specific column (by header name) for specific rows (by 'r' column).
 *
 * @param filename       Path to the input text file.
 * @param column_name    Header name of the desired column (e.g., "xi10").
 * @param rows           Array of row identifiers to fetch (values from the 'r' column).
 * @param nrows          Number of entries in 'rows'.
 * @return               Pointer to a newly allocated array of doubles of length nrows, or NULL on error.
 *                       The i-th element corresponds to rows[i]. If a row isn't found, result[i] = NAN.
 *
 * Notes:
 *  - Caller owns the returned pointer and must free() it.
 *  - File is expected to be whitespace-delimited; header line may start with '#'.
 *  - Lines beginning with '#' (after leading whitespace) are treated as comments and skipped,
 *    except the first header line which is parsed even if it starts with '#'.
 */
double* read_values_by_rows(const char *filename,
                            const char *column_name,
                            const double *rows,
                            size_t nrows,
                            double *out)
{
    const char *delims = " \t\r\n";
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        perror("fopen");
        return NULL;
    }

    char line[8192];

    // --- 1) Read header line (first non-empty line). If it starts with '#', strip it and parse.
    int col_idx = -1;        // index of target column within the row (r=0, xi1=1, ...)
    int parsed_header = 0;

    while (fgets(line, sizeof(line), fp)) {
        // Trim leading whitespace
        char *p = line;
        while (*p == ' ' || *p == '\t') p++;
        // Skip empty lines
        if (*p == '\0' || *p == '\n' || *p == '\r') continue;

        // If line begins with '#', this can be the header line; strip one leading '#'
        if (*p == '#') p++;

        // Make a copy for tokenization
        char header_copy[8192];
        strncpy(header_copy, p, sizeof(header_copy));
        header_copy[sizeof(header_copy) - 1] = '\0';

        // Tokenize header to discover column indices
        char *tok = strtok(header_copy, delims);
        int idx = 0;
        while (tok) {
            if (strcmp(tok, column_name) == 0) {
                col_idx = idx;
            }
            tok = strtok(NULL, delims);
            idx++;
        }

        if (col_idx < 0) {
            fprintf(stderr, "Error: Column '%s' not found in header.\n", column_name);
            fclose(fp);
            return NULL;
        }

        parsed_header = 1;
        break;  // header parsed
    }

    if (!parsed_header) {
        fprintf(stderr, "Error: Could not read header line.\n");
        fclose(fp);
        return NULL;
    }

    // Prepare output; initialize to NAN to flag missing rows.
    // double *out = (double *)malloc(nrows * sizeof(double));
    // if (!out) {
    //     perror("malloc");
    //     fclose(fp);
    //     return NULL;
    // }
    for (size_t i = 0; i < nrows; i++) {
        out[i] = NAN;
    }

    // --- 2) Read data rows
    while (fgets(line, sizeof(line), fp)) {
        // Trim leading whitespace
        char *p = line;
        while (*p == ' ' || *p == '\t') p++;

        // Skip empty or comment lines
        if (*p == '\0' || *p == '\n' || *p == '\r' || *p == '#') continue;

        // Tokenize the line
        char row_copy[8192];
        strncpy(row_copy, p, sizeof(row_copy));
        row_copy[sizeof(row_copy) - 1] = '\0';

        int idx = 0;
        int row_id = 0;
        double value = 0.0;
        int have_row = 0, have_value = 0;

        char *tok = strtok(row_copy, delims);
        while (tok) {
            if (idx == 0) {
                // First column is 'r'
                char *endp = NULL;
                row_id = (int)strtol(tok, &endp, 10);
                if (endp == tok) break; // not a number; malformed row
                have_row = 1;
            }
            if (idx == col_idx) {
                char *endp = NULL;
                value = strtod(tok, &endp);
                if (endp == tok) {
                    // malformed value; treat as missing
                    have_value = 0;
                } else {
                    have_value = 1;
                }
            }
            idx++;
            tok = strtok(NULL, delims);
        }

        if (!have_row || !have_value) continue;

        // Assign value to every matching request index (keeps input order, supports duplicates)
        for (size_t i = 0; i < nrows; i++) {
            if (llround(rows[i]) == row_id) {
                out[i] = value;
            }
        }
    }

    fclose(fp);
}