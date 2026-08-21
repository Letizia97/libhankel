// interp_cubic.cpp
//
// C++ implementation behind the C API in interp_cubic.h, same pattern as
// boost_bessel.cpp: Boost and the C++ types stay out of the header, and the
// API functions are extern "C" so they link under their plain names.

// WORKAROUND (Boost 1.74): pchip.hpp's constructor calls plain isnan() (lines
// 28, 59) but only says `using std::isnan` inside push_back() (line 86), so it
// will not compile. The fix must come *before* the include -- the name is
// resolved as the template is read. Delete both lines once Boost is past 1.74.
#include <cmath>
using std::isnan;

#include <boost/math/interpolators/pchip.hpp>

#include <exception>
#include <vector>

#include "interp_cubic.h"

// The definition the header withheld -- what makes cubic_interp_t opaque to C.
// The tag must match the header's typedef. Held by value, so destroying the
// struct destroys the spline. Outside extern "C" by convention: a linkage
// specification only applies to function and variable names, so it would do
// nothing here.
struct cubic_interp {
    boost::math::interpolators::pchip<std::vector<double>> spline;
};

// C cannot catch a C++ exception: one escaping any function below kills the
// process on the spot (from Python, the interpreter vanishes with no
// traceback). So each catches its own errors and reports failure by return
// value.
extern "C" {

cubic_interp_t *cubic_interp_create(const double *x, const double *y, size_t n) {
    // Outside the try on purpose: reading through a null pointer is a segfault,
    // not an exception, so the catch below could not turn it into a NULL.
    if (x == NULL || y == NULL)
        return NULL;

    try {
        // Shorthand for boost::math::interpolators::pchip.
        using boost::math::interpolators::pchip;

        // Boost holds its input for the spline's lifetime, but x and y are
        // borrowed and die with the caller's frame, so copy into owning
        // containers.
        std::vector<double> xv(x, x + n);
        std::vector<double> yv(y, y + n);

        // The angle brackets take the container, not double -- Boost reads the
        // number type off the vector. std::move is required:
        // the constructor's `&&` parameters only accept values it may take 
        // ownership of, and a named variable does not qualify. xv and
        // yv are emptied here. If the constructor throws, the `new` allocation
        // is released, so the catch below does not leak.
        return new cubic_interp{pchip<std::vector<double>>(std::move(xv), std::move(yv))};
    } catch (const std::exception &) {
        // domain_error for n < 4 or x not strictly increasing; bad_alloc from
        // `new`. NULL is the C convention for "construction failed".
        return NULL;
    }
}

double cubic_interp_eval(const cubic_interp_t *h, double xi) {
    // A null deref is a segfault, not something the catch below could see.
    if (h == NULL)
        return NAN;

    try {
        return h->spline(xi);
    } catch (const std::exception &) {
        // Boost refuses to extrapolate: xi outside [x[0], x[n-1]] throws.
        // NaN matches boost_bessel.cpp's convention.
        return NAN;
    }
}

void cubic_interp_destroy(cubic_interp_t *h) {
    // No try/catch needed: destructors cannot throw. delete on NULL is a no-op.
    delete h;
}
}
