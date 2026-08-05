#ifndef VALIDATE_X_H
#define VALIDATE_X_H

#include <stddef.h>

/**
 * @brief Rejects a transform variable that is not finite and greater than zero.
 *
 * Every strategy places its nodes at some constant divided by x, so a value
 * that is zero, negative, NaN or infinite would be reported as a success
 * carrying Inf, NaN, DBL_MAX or - for x < 0 - a plausible-looking but wrong
 * number.  Shared by all strategies so the condition, the message and the
 * status code cannot drift apart between them.
 *
 * @param x  value to check
 *
 * @return 0 if @p x is usable, or -12 after writing a message to stderr.
 */
int validate_x(double x);

/**
 * @brief Applies @ref validate_x to every entry of an array.
 *
 * Checked up front rather than per point, so an invalid input costs nothing:
 * the strategies that take a whole array precompute their node tables before
 * the loop over points, and there is no reason to build them for a call that
 * cannot succeed.
 *
 * @param x      array to check
 * @param len_x  number of entries in @p x
 *
 * @return 0 if every entry is usable, or -12 for the first one that is not.
 */
int validate_x_array(const double *x, size_t len_x);

#endif // VALIDATE_X_H
