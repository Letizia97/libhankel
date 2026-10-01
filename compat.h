#ifndef LIBHANKEL_COMPAT_H
#define LIBHANKEL_COMPAT_H

/* M_PI and M_PI_2 are POSIX extensions, not C standard. Define them here if
 * the platform's <math.h> did not (e.g. MSVC without _USE_MATH_DEFINES). */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

#endif /* LIBHANKEL_COMPAT_H */
