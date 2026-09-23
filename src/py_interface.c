
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
static double* python_sequence_to_c_array(PyObject *seq_obj, Py_ssize_t *out_len) {
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

static PyObject *py_hankel_transform(PyObject *self, PyObject *args) {
    int nu;
    PyObject *f_obj, *x_obj, *params_obj, *strategy_param_obj;
    const char *strategy_name;

    // Parse: nu, form_factor, x_points, form_factor_params,
    // strategy_name, strategy_params_dict
    if (!PyArg_ParseTuple(args, "iOOOsO", &nu, &f_obj, &x_obj, &params_obj, &strategy_name,
                          &strategy_param_obj)) {
        return NULL;
    }

    // Initialize all pointers to NULL for cleanup tracking.
    double *x = NULL;
    double *f_params = NULL;
    double *output = NULL;
    py_f_ctx *f_ctx = NULL;
    form_factor_f f_ptr = NULL;

    // ---------------------------
    // Convert x → C array
    // ---------------------------
    Py_ssize_t len_x;
    x = python_sequence_to_c_array(x_obj, &len_x);
    if (!x) {
        goto cleanup;
    }

    // ---------------------------
    // Convert params → C array
    // ---------------------------
    Py_ssize_t n_params;
    f_params = python_sequence_to_c_array(params_obj, &n_params);
    if (!f_params) {
        goto cleanup;
    }

    // ---------------------------
    // strategy_params struct
    // ---------------------------
    if (!PyDict_Check(strategy_param_obj)) {
        PyErr_SetString(PyExc_TypeError, "strategy_params must be dict");
        goto cleanup;
    }

    strategy_params sp;
    PyObject *n_eval_obj = PyDict_GetItemString(strategy_param_obj, "n_eval");
    PyObject *eps_rel_obj = PyDict_GetItemString(strategy_param_obj, "eps_rel");
    PyObject *f_max_obj = PyDict_GetItemString(strategy_param_obj, "f_max");

    if (!n_eval_obj) {
        sp.n_eval = 0;
    } else {
        sp.n_eval = PyFloat_AsDouble(n_eval_obj);
    }

    if (!eps_rel_obj) {
        sp.eps_rel = 0.0;
    } else {
        sp.eps_rel = PyFloat_AsDouble(eps_rel_obj);
    }

    if (!f_max_obj) {
        sp.f_max = 0.0;
    } else {
        sp.f_max = PyFloat_AsDouble(f_max_obj);
    }

    if (PyErr_Occurred()) {
        goto cleanup;
    }

    // ---------------------------
    // Output allocation
    // ---------------------------
    output = malloc(len_x * sizeof(double));
    if (!output) {
        PyErr_SetString(PyExc_MemoryError, "alloc output failed");
        goto cleanup;
    }

    // ---------------------------
    // Form factor dispatch
    // ---------------------------
    // f_obj can be:
    //   1. A Python callable (user function) → python_form_factor wrapper
    //   2. A string naming a built-in (e.g. 'gdab') → lookup by name
    //   3. A TabulatedFormFactor instance → [TO BE ADDED]
    f_ctx = malloc(sizeof(py_f_ctx));
    if (!f_ctx) {
        PyErr_SetString(PyExc_MemoryError, "f_ctx alloc failed");
        goto cleanup;
    }

    f_ctx->params = f_params;
    f_ctx->n_params = n_params;
    f_ctx->callable = NULL;

    // Case 1: Python callable (user-defined function)
    if (PyCallable_Check(f_obj)) {
        f_ptr = python_form_factor;
        f_ctx->callable = f_obj;
        Py_INCREF(f_obj);

    // Case 2: Built-in form factor by string name
    } else if (PyUnicode_Check(f_obj)) {
        const char *name = PyUnicode_AsUTF8(f_obj);
        f_ptr = get_form_factor_by_name(name);

        if (!f_ptr) {
            PyErr_SetString(PyExc_ValueError, "Unknown function name");
            goto cleanup;
        }

    // Case 3: [Will add TabulatedFormFactor instance detection here]
    // Need to check whether data as been passed in
    // then call interp on the data and then get 

    } else if (PyDict_Check(f_obj)) {
        // Extract Python objects from dict
        PyObject *q_obj = PyDict_GetItemString(f_obj, "q");
        PyObject *f_obj_data = PyDict_GetItemString(f_obj, "f");
        PyObject *interp_type_obj = PyDict_GetItemString(f_obj, "interp_type");
        PyObject *tail_obj = PyDict_GetItemString(f_obj, "tail");

        PyObject *exponent_obj = PyDict_GetItemString(f_obj, "exponent");
        double exponent = 0.0;
        if (exponent_obj) {
            exponent = PyFloat_AsDouble(exponent_obj);
        }

        double *q_array = NULL;
        double *f_array = NULL;

        // Convert to C arrays
        Py_ssize_t len_q;
        double *q_array = python_sequence_to_c_array(q_obj, &len_q);
        if (!q_array) {
            goto cleanup;
        }

        Py_ssize_t len_f;
        double *f_array = python_sequence_to_c_array(f_obj_data, &len_f);
        if (!f_array) {
            free(q_array);
            goto cleanup;
        }

        // Map the string interp_type to enum values
        const char *interp_type_str = PyUnicode_AsUTF8(interp_type_obj);
        tabulated_interp_type_t interp_type;
        if (strcmp(interp_type_str, "linear") == 0) {
            interp_type = TABULATED_INTERP_LINEAR;
        } else if (strcmp(interp_type_str, "cubic") == 0) {
            interp_type = TABULATED_INTERP_CUBIC;
        } else if (strcmp(interp_type_str, "loglinear") == 0) {
            interp_type = TABULATED_INTERP_LOGLINEAR;
        } else {
            PyErr_SetString(PyExc_ValueError, "unknown interp_type");
            goto cleanup;
        } 

        // Map the string tail to enum values
        const char *tail_str = PyUnicode_AsUTF8(tail_obj);
        tabulated_tail_t tail;
        if (strcmp(tail_str, "power_law") == 0) {
            tail = TABULATED_TAIL_POWER_LAW;
        } else if (strcmp(tail_str, "zero") == 0) {
            tail = TABULATED_TAIL_ZERO;
        } else {
            PyErr_SetString(PyExc_ValueError, "unknown tail");
            goto cleanup;
        } 

        tabulated_ff_t *tff_handle = NULL;
        int status = tabulated_ff_create(q_array, f_array, len_q,
                                        interp_type, tail, exponent,
                                        &tff_handle);
        if (status != 0) {
            PyErr_SetString(
                PyExc_ValueError, "Failed to create tabulated form factor"
            );
            free(q_array);
            free(f_array);
            goto cleanup;
        }
        
        // points to the interpolation evaluator
        f_ptr = tabulated_ff_eval;

        // opaque struct with spline data
        f_ctx = (void *)tff_handle;

    } else {
        PyErr_SetString(PyExc_TypeError, "f must be callable, string, or TabulatedFormFactor");
        goto cleanup;
    }

    int status_code = hankel_transform(nu, f_ptr, x, len_x, f_ctx, output, strategy_name, sp);

    // Error handling
    switch (status_code) {
    case 0:
        break;

    case -1:
        PyErr_SetString(PyExc_ValueError,
                        "nu needs to be 0 or 1 in order to use the selected strategy");
        goto cleanup;

    case -2:
        /* Unreachable from Python: strategies are selected by name here, and
         * hankel_transform() only ever passes a valid index to the filters. */
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

    // Build output list on success
    PyObject *out_list = PyList_New(len_x);
    for (Py_ssize_t i = 0; i < len_x; i++) {
        PyList_SetItem(out_list, i, PyFloat_FromDouble(output[i]));
    }

    // Clean up and return success
    free(x);
    free(output);
    if (tff_handle) {
        tabulated_ff_destroy(tff_handle);
    }
    free(f_ctx->params);
    if (f_ctx->callable) {
        Py_DECREF(f_ctx->callable);
    }
    free(f_ctx);
    return out_list;

cleanup:
    // On error, free all allocated resources
    free(x);
    free(output);
    if (tff_handle) {
        tabulated_ff_destroy(tff_handle);
    }
    if (f_ctx) {
        free(f_ctx->params);
        if (f_ctx->callable) {
            Py_DECREF(f_ctx->callable);
        }
        free(f_ctx);
    } else {
        // f_ctx was never allocated, so f_params wasn't transferred to it
        free(f_params);
    }
    return NULL;
}

// docstring for hankel_tranform python api
static char hankel_t_doc[] =
    "Compute the Hankel transform.\n"
    "\n"
    ":param nu:              The order of bessel function, must be 0 or 1.\n"
    ":type nu:               int \n"
    ":param f:               Either a function to hankel-transform or a string "
    "naming a built-in function, "
    "i.e. either of 'gdab', 'broad_peak', 'sphere'.\n"
    ":type f:                callable or str\n"
    ":param x_arr:           The points at which to evaluate the Hankel "
    "transform of the function f.\n"
    ":type x_arr:            numpy.ndarray of float64\n"
    ":param f_params:        Input parameters needed by the function f "
    "(ordered).\n"
    ":type f_params:         numpy.ndarray of float64\n"
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
    "function with either a "
    "builtin form factor or a custom input function."
    "\n";

static PyMethodDef Methods[] = {
    {"hankel_transform", py_hankel_transform, METH_VARARGS, hankel_t_doc}, {NULL, NULL, 0, NULL}};

static struct PyModuleDef module = {PyModuleDef_HEAD_INIT, "libhankel", " ", -1, Methods};

PyMODINIT_FUNC PyInit_libhankel(void) { return PyModule_Create(&module); }
