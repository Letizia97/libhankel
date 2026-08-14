// interp_cubic.cpp
//
// C++ implementation behind the C API declared in interp_cubic.h.
// Same pattern as boost_bessel.cpp: the Boost headers and every C++ type
// stay in this file, and only unmangled C symbols escape.

// ---------------------------------------------------------------------------
// WORKAROUND -- the ordering below is load-bearing. Do not tidy it.
//
// Boost 1.74's pchip.hpp calls isnan() unqualified in the constructor
// (lines 28 and 59) but only declares `using std::isnan` inside push_back()
// (line 86). The argument is a double, a built-in type with no associated
// namespace, so argument-dependent lookup finds nothing, and <cmath> puts
// isnan in std:: rather than the global namespace. Result: "'isnan' was not
// declared in this scope".
//
// Hoisting std::isnan into the global namespace fixes it, but the using
// declaration MUST come before the pchip.hpp include: the call is resolved
// by ordinary lookup at the point the template is *defined*, i.e. while
// pchip.hpp is being parsed. Placing it after has no effect.
//
// Delete both lines once Boost is upgraded past 1.74.
#include <cmath>
using std::isnan;
// ---------------------------------------------------------------------------

#include <boost/math/interpolators/pchip.hpp>

#include <exception>
#include <vector>

#include "interp_cubic.h"

// The definition the header deliberately withheld. Keeping it here is what
// makes cubic_interp_t opaque to C.
//
// The tag name must match the header's `typedef struct cubic_interp` exactly,
// or the type stays incomplete and every use fails.
//
// The spline is held by value, not behind a std::unique_ptr. The struct
// itself is what gets heap-allocated, and pchip cleans up after itself when
// the struct is destroyed -- one level of indirection instead of two.
//
// Outside extern "C" on purpose: that block controls how *function* names are
// mangled for the linker, and a struct definition has no linkage to specify.
struct cubic_interp {
    boost::math::interpolators::pchip<std::vector<double>> spline;
};

// Everything below is a C boundary. No exception may escape any of these
// functions: C frames carry no unwinding information, so an exception
// propagating out of one reaches no handler and the runtime calls
// std::terminate() -- an immediate SIGABRT with no chance to recover. From
// Python that kills the interpreter outright, with no traceback.
extern "C" {

cubic_interp_t *cubic_interp_create(const double *x, const double *y, size_t n) {
    try {
        // The template parameter is the *container*, not the scalar. Boost
        // derives the scalar itself via RandomAccessContainer::value_type
        // (pchip.hpp:17), so pchip<double> would expand to double::value_type
        // and fail with a wall of template errors.
        using boost::math::interpolators::pchip;

        // x and y are bare pointers: they do not own their memory, do not
        // know their length, and will not outlive the caller's frame. Boost
        // stores its input for the lifetime of the spline, so it needs real
        // containers. These constructors take a half-open iterator range
        // [first, last), and a raw pointer is a valid iterator -- so x + n
        // marks one past the end and exactly n elements are copied.
        std::vector<double> xv(x, x + n);
        std::vector<double> yv(y, y + n);

        // std::move is mandatory here, not an optimisation. Boost's
        // constructor takes RandomAccessContainer&& -- an rvalue reference,
        // which only binds to temporaries. xv is a named variable (an
        // lvalue), so passing it bare is a hard compile error, not a silent
        // copy. std::move does no work at runtime; it just recategorises the
        // expression as an rvalue, i.e. "I am finished with this, take its
        // buffer". xv and yv are empty afterwards -- do not reuse them.
        //
        // Note pchip<...>(args) is an *expression* producing a temporary.
        // pchip<...> spline(args) would be a declaration, and you cannot
        // return a declaration.
        //
        // If the pchip constructor throws part-way through, C++ releases the
        // memory `new` just allocated, so the catch below does not leak.
        return new cubic_interp{pchip<std::vector<double>>(std::move(xv), std::move(yv))};
    } catch (const std::exception &) {
        // Boost throws std::domain_error for n < 4 (pchip.hpp:23) and for x
        // not strictly increasing (cubic_hermite_detail.hpp:50); `new` throws
        // std::bad_alloc. NULL is the C convention for "construction failed"
        // -- note this function returns a pointer, so NAN would not even
        // compile here.
        return NULL;
    }
}

double cubic_interp_eval(const cubic_interp_t *h, double xi) {
    // Guards against a caller who ignored a NULL from create(). A null
    // dereference is a segfault, not an exception, so the catch below would
    // never see it.
    if (h == NULL)
        return NAN;

    try {
        // h is a pointer to the struct, not to the spline -- reach the member
        // through it. h->spline is shorthand for (*h).spline.
        return h->spline(xi);
    } catch (const std::exception &) {
        // Boost range-checks and refuses to extrapolate: xi outside
        // [x[0], x[n-1]] throws std::domain_error (cubic_hermite_detail.hpp:77).
        // NaN matches what the wrappers in boost_bessel.cpp return, and
        // callers detect it with isnan().
        return NAN;
    }
}

void cubic_interp_destroy(cubic_interp_t *h) {
    // No try/catch here, unlike the two above. delete runs the struct's
    // destructor, and destructors are implicitly noexcept since C++11 -- if
    // one did throw, std::terminate() fires immediately and no handler could
    // catch it. Releasing a shared_ptr and its vectors cannot fail anyway, so
    // a catch block would be dead code implying a risk that does not exist.
    //
    // delete on a null pointer is a guaranteed no-op, so destroy(NULL) is safe.
    delete h;
}
}
