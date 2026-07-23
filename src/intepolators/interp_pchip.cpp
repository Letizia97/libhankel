#include "interpolation.h"
#include <boost/math/interpolators/pchip.hpp>

struct pchip_interp {
    // ptr here needed so that when interp is destroyed the
    // pchip_interp object is destroyed too
    std::unique_ptr<boost::math::interpolators::pchip<double>> pchip_interp;
};

pchip_interp *pchip_create(const double *x, const double *y, size_t n) {
    // allocate the implementation object
    auto *impl = new pchip_interp;

    // creates vectors from raw C arrays
    auto xv = std::vector<double>(x, x + n);
    auto yv = std::vector<double>(y, y + n);

    // access the interp member inside struct
    // this line means (*impl).interp
    impl->interp =
        // creates a new PCHIP interpolator object
        std::make_unique<boost::math::interpolators::pchip<double>
                         // Without move Boost would copy the vectors
                         // while using move allows Boost to "steal"/transfer their content
                         // Only used here because we don't need xv and yv after
                         >(std::move(xv), std::move(yv));

    return impl;
}

double pchip_eval(pchip_interp *impl, double x) {
    // auto = "figure out the type for me"
    //  & = "make it a reference"
    auto &interp = *(impl->interp);
    return interp(x);
}

void pchip_destroy(pchip_interp *impl) { delete impl; }
