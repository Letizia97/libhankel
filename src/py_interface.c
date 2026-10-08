
#define PY_SSIZE_T_CLEAN
#include "form_factors.h"
#include "libhankel.h"
#include <Python.h>
#include <stdlib.h>
#include "tabulated_ff.h"

typedef struct {
    double *params;
    size_t n_params;
    PyObject *callable; // NULL for built-ins
} py_f_ctx;

double python_form_factor(double x, void *f_ctx) {
    // Callback wrapper: invoke a user-defined Python function from C code.
    // Called by hankel_transform() for each evaluation point.
    // f_ctx holds the Python callable and its parameter list.
    py_f_ctx *c = (py_f_ctx *)f_ctx;

    // Acquire GIL: Python object manipulation requires the lock.
    PyGILState_STATE gstate = PyGILState_Ensure();

    // Build args tuple: (x, [param1, param2, ...])
    PyObject *args = PyTuple_New(2);
    PyTuple_SetItem(args, 0, PyFloat_FromDouble(x));

    // Pack parameters into a list.
    PyObject *list = PyList_New(c->n_params);
    for (size_t i = 0; i < c->n_params; i++) {
        PyList_SetItem(list, i, PyFloat_FromDouble(c->params[i]));
    }
    PyTuple_SetItem(args, 1, list);

    // Call the Python function with (x, params).
    PyObject *result = PyObject_CallObject(c->callable, args);
    Py_DECREF(args);

    // Extract the float result, or 0 on error.
    double val = 0.0;
    if (result) {
        val = PyFloat_AsDouble(result);
        Py_DECREF(result);
    } else {
        PyErr_Print();
    }

    PyGILState_Release(gstate);
    return val;
}

// Convert a Python sequence of floats to a C array. Caller must free the result.
// On error, sets a Python exception and returns NULL.
static double *python_sequence_to_c_array(PyObject *seq_obj, Py_ssize_t *out_len) {
    if (!PySequence_Check(seq_obj)) {
        PyErr_SetString(PyExc_TypeError, "expected a sequence");
        return NULL;
    }

    Py_ssize_t len = PySequence_Size(seq_obj);
    if (len < 0) {
        return NULL;
    }

    double *arr = malloc(len * sizeof(double));
    if (!arr) {
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate array");
        return NULL;
    }

    for (Py_ssize_t i = 0; i < len; i++) {
        PyObject *item = PySequence_GetItem(seq_obj, i);
        if (!item) {
            free(arr);
            return NULL;
        }

        arr[i] = PyFloat_AsDouble(item);
        Py_DECREF(item);

        if (PyErr_Occurred()) {
            free(arr);
            return NULL;
        }
    }

    *out_len = len;
    return arr;
}

// Resolved form factor: owns all resources, freed by free_resolved_ff().
typedef struct {
    form_factor_f f_ptr;
    void *ctx;           // passed to f_ptr; points at py_ctx or tff
    py_f_ctx *py_ctx;    // non-NULL for callable / built-in paths
    tabulated_ff_t *tff; // non-NULL for tabulated dict path
    double *f_params;    // always owned here; points into py_ctx->params when py_ctx != NULL
} resolved_ff;

// Parse a Python form-factor argument into a resolved_ff.
// f_obj:      callable, str, or dict
// params_obj: sequence of floats (empty list for tabulated)
// Returns 0 on success, -1 with a Python exception set on failure.
// Always call free_resolved_ff() after a successful call.
static int resolve_form_factor(PyObject *f_obj, PyObject *params_obj, resolved_ff *out) {
    *out = (resolved_ff){0};

    Py_ssize_t n_params;
    out->f_params = python_sequence_to_c_array(params_obj, &n_params);
    if (!out->f_params) {
        return -1;
    }

    if (PyCallable_Check(f_obj)) {
        out->py_ctx = malloc(sizeof(py_f_ctx));
        if (!out->py_ctx) {
            PyErr_NoMemory();
            return -1;
        }
        out->py_ctx->params = out->f_params;
        out->py_ctx->n_params = n_params;
        out->py_ctx->callable = f_obj;
        Py_INCREF(f_obj);
        out->f_ptr = python_form_factor;
        out->ctx = out->py_ctx;

    } else if (PyUnicode_Check(f_obj)) {
        out->py_ctx = malloc(sizeof(py_f_ctx));
        if (!out->py_ctx) {
            PyErr_NoMemory();
            return -1;
        }
        out->py_ctx->params = out->f_params;
        out->py_ctx->n_params = n_params;
        out->py_ctx->callable = NULL;
        const char *name = PyUnicode_AsUTF8(f_obj);
        out->f_ptr = get_form_factor_by_name(name);
        if (!out->f_ptr) {
            PyErr_SetString(PyExc_ValueError, "Unknown function name");
            return -1;
        }
        out->ctx = out->py_ctx;

    } else if (PyDict_Check(f_obj)) {
        PyObject *q_obj          = PyDict_GetItemString(f_obj, "q");
        PyObject *f_obj_data     = PyDict_GetItemString(f_obj, "f");
        PyObject *interp_type_obj = PyDict_GetItemString(f_obj, "interp_type");
        PyObject *tail_obj       = PyDict_GetItemString(f_obj, "tail");

        if (!q_obj || !f_obj_data || !interp_type_obj || !tail_obj) {
            PyErr_SetString(PyExc_ValueError,
                "tabulated dict must contain 'q', 'f', 'interp_type', and 'tail'");
            return -1;
        }

        PyObject *exponent_obj = PyDict_GetItemString(f_obj, "exponent");
        double exponent = 0.0;
        if (exponent_obj) {
            exponent = PyFloat_AsDouble(exponent_obj);
        }

        Py_ssize_t len_q;
        double *q_array = python_sequence_to_c_array(q_obj, &len_q);
        if (!q_array) {
            return -1;
        }

        Py_ssize_t len_f;
        double *f_array = python_sequence_to_c_array(f_obj_data, &len_f);
        if (!f_array) {
            free(q_array);
            return -1;
        }

        const char *interp_type_str = PyUnicode_AsUTF8(interp_type_obj);
        tabulated_interp_type_t interp_type;
        if (strcmp(interp_type_str, "linear") == 0) {
            interp_type = TABULATED_INTERP_LINEAR;
        } else if (strcmp(interp_type_str, "cubic") == 0) {
            interp_type = TABULATED_INTERP_CUBIC;
        } else if (strcmp(interp_type_str, "loglinear") == 0) {
            interp_type = TABULATED_INTERP_LOGLINEAR;
        } else {
            free(q_array);
            free(f_array);
            PyErr_SetString(PyExc_ValueError,
                "interp_type must be 'linear', 'cubic', or 'loglinear'");
            return -1;
        }

        const char *tail_str = PyUnicode_AsUTF8(tail_obj);
        tabulated_tail_t tail;
        if (strcmp(tail_str, "power_law") == 0) {
            tail = TABULATED_TAIL_POWER_LAW;
        } else if (strcmp(tail_str, "zero") == 0) {
            tail = TABULATED_TAIL_ZERO;
        } else {
            free(q_array);
            free(f_array);
            PyErr_SetString(PyExc_ValueError, "tail must be 'power_law' or 'zero'");
            return -1;
        }

        int status = tabulated_ff_create(q_array, f_array, len_q,
                                         interp_type, tail, exponent, &out->tff);
        free(q_array);
        free(f_array);
        if (status != 0) {
            PyErr_SetString(PyExc_ValueError,
                "Failed to create tabulated form factor (check q/f ranges and types)");
            return -1;
        }

        out->f_ptr = tabulated_ff_eval;
        out->ctx   = out->tff;

    } else {
        PyErr_SetString(PyExc_TypeError, "f must be callable, string, or dict");
        return -1;
    }

    return 0;
}

static void free_resolved_ff(resolved_ff *ff) {
    if (ff->tff) {
        tabulated_ff_destroy(ff->tff);
    }
    if (ff->py_ctx) {
        if (ff->py_ctx->callable) {
            Py_DECREF(ff->py_ctx->callable);
        }
        free(ff->py_ctx);
    }
    free(ff->f_params);
}

static PyObject *py_hankel_transform(PyObject *self, PyObject *args) {
    int nu;
    PyObject *f_obj, *x_obj, *params_obj, *strategy_param_obj;
    const char *strategy_name;

    if (!PyArg_ParseTuple(args, "iOOOsO", &nu, &f_obj, &x_obj, &params_obj, &strategy_name,
                          &strategy_param_obj)) {
        return NULL;
    }

    double *x = NULL;
    double *output = NULL;
    resolved_ff ff = {0};
    PyObject *out_list = NULL;

    Py_ssize_t len_x;
    x = python_sequence_to_c_array(x_obj, &len_x);
    if (!x) {
        goto cleanup;
    }

    if (!PyDict_Check(strategy_param_obj)) {
        PyErr_SetString(PyExc_TypeError, "strategy_params must be dict");
        goto cleanup;
    }

    strategy_params sp;
    PyObject *n_eval_obj = PyDict_GetItemString(strategy_param_obj, "n_eval");
    PyObject *eps_rel_obj = PyDict_GetItemString(strategy_param_obj, "eps_rel");
    PyObject *f_max_obj = PyDict_GetItemString(strategy_param_obj, "f_max");

    sp.n_eval  = n_eval_obj  ? (int)PyFloat_AsDouble(n_eval_obj)  : 0;
    sp.eps_rel = eps_rel_obj ? PyFloat_AsDouble(eps_rel_obj)       : 0.0;
    sp.f_max   = f_max_obj   ? PyFloat_AsDouble(f_max_obj)         : 0.0;

    if (PyErr_Occurred()) {
        goto cleanup;
    }

    output = malloc(len_x * sizeof(double));
    if (!output) {
        PyErr_SetString(PyExc_MemoryError, "alloc output failed");
        goto cleanup;
    }

    if (resolve_form_factor(f_obj, params_obj, &ff) < 0) {
        goto cleanup;
    }

    int status_code = hankel_transform(nu, ff.f_ptr, x, len_x, ff.ctx, output, strategy_name, sp);

    switch (status_code) {
    case 0:
        break;
    case -1:
        PyErr_SetString(PyExc_ValueError,
                        "nu needs to be 0 or 1 in order to use the selected strategy");
        goto cleanup;
    case -2:
        PyErr_SetString(PyExc_RuntimeError, "Internal error: invalid DHT filter index");
        goto cleanup;
    case -3:
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate internal variables");
        goto cleanup;
    case -4:
        PyErr_SetString(PyExc_RuntimeError, "Failed to converge");
        goto cleanup;
    case -5:
        PyErr_SetString(PyExc_ZeroDivisionError, "Internal error: division by zero");
        goto cleanup;
    case -6:
        PyErr_SetString(PyExc_ValueError,
                        "Internal error: wrong nzeros in function bessel_j_zero (must be >= 1)");
        goto cleanup;
    case -7:
        PyErr_SetString(PyExc_ValueError,
                        "Internal error: wrong n of iterations in pade sum (must be >= 1)");
        goto cleanup;
    case -8:
        PyErr_SetString(PyExc_ValueError, "Error: n_eval must be provided and cannot be zero");
        goto cleanup;
    case -9:
        PyErr_SetString(PyExc_ValueError, "Error: eps_rel must be provided and cannot be zero");
        goto cleanup;
    case -10:
        PyErr_SetString(PyExc_ValueError, "Error: f_max must be provided and cannot be zero");
        goto cleanup;
    case -11:
        PyErr_SetString(PyExc_ValueError,
                        "Error: invalid strategy name, must be one of : " LIBHANKEL_ALL_STRATEGIES
                        ".");
        goto cleanup;
    case -12:
        PyErr_SetString(PyExc_ValueError, "Error: x must be finite and greater than zero");
        goto cleanup;
    default:
        PyErr_SetString(PyExc_RuntimeError, "unknown error");
        goto cleanup;
    }

    out_list = PyList_New(len_x);
    for (Py_ssize_t i = 0; i < len_x; i++) {
        PyList_SetItem(out_list, i, PyFloat_FromDouble(output[i]));
    }

cleanup:
    free(x);
    free(output);
    free_resolved_ff(&ff);
    return out_list;
}

static PyObject *py_hankel_G0(PyObject *self, PyObject *args) {
    PyObject *f_obj, *params_obj;
    double epsrel;

    if (!PyArg_ParseTuple(args, "OOd", &f_obj, &params_obj, &epsrel)) {
        return NULL;
    }

    resolved_ff ff = {0};
    PyObject *result = NULL;

    if (resolve_form_factor(f_obj, params_obj, &ff) < 0) {
        goto cleanup;
    }

    result = PyFloat_FromDouble(hankel_G0(ff.f_ptr, ff.ctx, epsrel));

cleanup:
    free_resolved_ff(&ff);
    return result;
}

static char hankel_G0_doc[] =
    "Compute G(0), the Hankel transform evaluated at r = 0.\n"
    "\n"
    "Computes the integral of q * f(q) from 0 to infinity.\n"
    "\n"
    ":param f:        Either a callable, a string naming a built-in form factor\n"
    "                ('gdab', 'broad_peak', 'sphere'), or a dict for tabulated\n"
    "                form factor with keys 'q', 'f', 'interp_type', 'tail',\n"
    "                and optional 'exponent'.\n"
    ":type f:         callable, str, or dict\n"
    ":param f_params: Parameters for built-in or callable form factors.\n"
    "                Pass an empty list for tabulated form factors.\n"
    ":type f_params:  list or numpy.ndarray\n"
    ":param epsrel:   Relative error tolerance (e.g. 1e-6).\n"
    ":type epsrel:    float\n"
    ":returns:        The value of G(0).\n"
    ":rtype:          float\n";

static char hankel_t_doc[] =
    "Compute the Hankel transform.\n"
    "\n"
    ":param nu:              The order of bessel function, must be 0 or 1.\n"
    ":type nu:               int \n"
    ":param f:               Either a callable, a string naming a built-in "
    "form factor ('gdab', 'broad_peak', 'sphere'), or a dict for tabulated "
    "form factor with keys 'q', 'f', 'interp_type', 'tail', and optional 'exponent'.\n"
    ":type f:                callable, str, or dict\n"
    ":param x_arr:           The points at which to evaluate the Hankel "
    "transform of the function f.\n"
    ":type x_arr:            numpy.ndarray of float64\n"
    ":param f_params:        Input parameters needed by built-in or callable form factors "
    "(ordered). Not needed for tabulated form factors, pass an empty list.\n"
    ":type f_params:         numpy.ndarray of float64 or list\n"
    ":param strategy_name:   The name of Hankel strategy to use. "
    "Refer to the table in :ref:`strategy-parameters` for a list of possible "
    "strategies.\n"
    ":type strategy_name:    str \n"
    ":param strategy_params: The parameters needed by the chosen strategy. "
    "See :ref:`strategy-selection` or :ref:`c-api` (strategy_params) for "
    "details.\n"
    ":type strategy_params:  dict[str, float | int]\n"
    ":returns:               The hankel transform.\n"
    ":rtype:                 numpy.ndarray of float64\n"
    "Please refer to :ref:`python-examples` for examples on how to use this "
    "function with callable, built-in, or tabulated form factors."
    "\n";

static PyMethodDef Methods[] = {
    {"hankel_transform", py_hankel_transform, METH_VARARGS, hankel_t_doc},
    {"hankel_G0", py_hankel_G0, METH_VARARGS, hankel_G0_doc},
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef module = {PyModuleDef_HEAD_INIT, "libhankel", " ", -1, Methods};

PyMODINIT_FUNC PyInit_libhankel(void) { return PyModule_Create(&module); }
