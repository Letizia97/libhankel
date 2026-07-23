#include "interpolation.h"

struct interp {
    interp_method_t method;
    void *impl;
};

interp_t *interp_create(...) {
    switch (method) {
    case INTERP_LINEAR:
        ... break;

    case INTERP_PCHIP: {
        interp->impl = pchip_create(x, y, n);
        break;
    }

    case INTERP_MAKIMA:
        ... break;
    }
}

double interp_eval(...) {
    switch (interp->method) { ... }
}