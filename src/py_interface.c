
#define PY_SSIZE_T_CLEAN
#include "form_factors.h"
#include "libhankel.h"
#include <Python.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    double *params;
    size_t n_params;
    PyObject *callable; // NULL for built-ins
} py_f_ctx;

/* Node tables for the fixed-abscissa digital filters, defined in hankel_DHT.c.
 * Column 0 is the abscissa, columns 1 and 2 the J0 and J1 weights. */
extern const double KK51Hankel[51][3];
extern const double KK101Hankel[101][3];
extern const double KK201Hankel[201][3];
extern const double WA801Hankel[801][3];

/* Look up the node table a fixed-abscissa filter uses, or NULL if the named
 * strategy chooses its abscissae adaptively and so has no table. */
static const double *dht_node_table(const char *strategy_name, size_t *n_nodes) {
    if (strcmp(strategy_name, "DHT_Key_51") == 0) {
        *n_nodes = 51;
        return &KK51Hankel[0][0];
    }
    if (strcmp(strategy_name, "DHT_Key_101") == 0) {
        *n_nodes = 101;
        return &KK101Hankel[0][0];
    }
    if (strcmp(strategy_name, "DHT_Key_201") == 0) {
        *n_nodes = 201;
        return &KK201Hankel[0][0];
    }
    if (strcmp(strategy_name, "DHT_Anderson_801") == 0) {
        *n_nodes = 801;
        return &WA801Hankel[0][0];
    }
    return NULL;
}

double python_form_factor(double x, void *f_ctx) {
    py_f_ctx *c = (py_f_ctx *)f_ctx;

    PyGILState_STATE gstate = PyGILState_Ensure();

    PyObject *args = PyTuple_New(2);
    PyTuple_SetItem(args, 0, PyFloat_FromDouble(x));

    PyObject *list = PyList_New(c->n_params);

    for (size_t i = 0; i < c->n_params; i++) {
        PyList_SetItem(list, i, PyFloat_FromDouble(c->params[i]));
    }

    PyTuple_SetItem(args, 1, list);

    PyObject *result = PyObject_CallObject(c->callable, args);
    Py_DECREF(args);

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

static PyObject *py_hankel_transform(PyObject *self, PyObject *args) {
    int nu;
    PyObject *f_obj, *x_obj, *params_obj, *strategy_param_obj;
    const char *strategy_name;

    if (!PyArg_ParseTuple(args, "iOOOsO", &nu, &f_obj, &x_obj, &params_obj, &strategy_name,
                          &strategy_param_obj)) {
        return NULL;
    }

    // ---------------------------
    // Convert x → C array
    // ---------------------------
    if (!PySequence_Check(x_obj)) {
        PyErr_SetString(PyExc_TypeError, "x must be a sequence");
        {
            return NULL;
        }
    }

    Py_ssize_t len_x = PySequence_Size(x_obj);
    if (len_x < 0) {
        return NULL;
    }

    double *x = malloc(len_x * sizeof(double));
    if (!x) {
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate x");
        return NULL;
    }

    for (Py_ssize_t i = 0; i < len_x; i++) {
        PyObject *item = PySequence_GetItem(x_obj, i);
        if (!item) {
            free(x);
            return NULL;
        }

        x[i] = PyFloat_AsDouble(item);
        Py_DECREF(item);

        if (PyErr_Occurred()) {
            free(x);
            return NULL;
        }
    }

    // ---------------------------
    // Convert params → C array
    // ---------------------------
    if (!PySequence_Check(params_obj)) {
        free(x);
        PyErr_SetString(PyExc_TypeError, "params must be a sequence");
        return NULL;
    }

    Py_ssize_t n_params = PySequence_Size(params_obj);
    double *f_params = malloc(n_params * sizeof(double));

    for (Py_ssize_t i = 0; i < n_params; i++) {
        PyObject *item = PySequence_GetItem(params_obj, i);
        if (!item) {
            free(x);
            return NULL;
        }

        f_params[i] = PyFloat_AsDouble(item);
        Py_DECREF(item);

        if (PyErr_Occurred()) {
            free(x);
            return NULL;
        }
    }

    // ---------------------------
    // strategy_params struct
    // ---------------------------
    if (!PyDict_Check(strategy_param_obj)) {
        free(x);
        PyErr_SetString(PyExc_TypeError, "strategy_params must be dict");
        return NULL;
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
        free(x);
        return NULL;
    }

    // ---------------------------
    // Output allocation
    // ---------------------------
    double *output = malloc(len_x * sizeof(double));
    if (!output) {
        free(x);
        PyErr_SetString(PyExc_MemoryError, "alloc output failed");
        return NULL;
    }

    // ---------------------------
    // Form factor function
    // ---------------------------
    form_factor_f f_ptr = NULL;
    py_f_ctx *f_ctx = malloc(sizeof(py_f_ctx));
    if (!f_ctx) {
        free(x);
        free(f_params);
        free(output);
        PyErr_SetString(PyExc_MemoryError, "f_ctx alloc failed");
        return NULL;
    }

    f_ctx->params = f_params;
    f_ctx->n_params = n_params;
    f_ctx->callable = NULL;

    /* check for Python callable */
    if (PyCallable_Check(f_obj)) {
        f_ptr = python_form_factor;
        f_ctx->callable = f_obj;
        Py_INCREF(f_obj);

        /* or whether the user wants to use built-in form factor */
    } else if (PyUnicode_Check(f_obj)) {
        const char *name = PyUnicode_AsUTF8(f_obj);
        f_ptr = get_form_factor_by_name(name);

        if (!f_ptr) {
            free(x);
            free(f_params);
            free(output);
            PyErr_SetString(PyExc_ValueError, "Unknown function name");
            return NULL;
        }

    } else {
        free(x);
        free(f_params);
        free(output);
        PyErr_SetString(PyExc_TypeError, "f must be callable or string");
        return NULL;
    }

    int status_code = hankel_transform(nu, f_ptr, x, len_x, f_ctx, output, strategy_name, sp);

    // Error handling
    switch (status_code) {
    case 0:
        break;

    case -1:
        PyErr_SetString(PyExc_ValueError,
                        "nu needs to be 0 or 1 in order to use the selected strategy");
        return NULL;

    case -2:
        /* Unreachable from Python: strategies are selected by name here, and
         * hankel_transform() only ever passes a valid index to the filters. */
        PyErr_SetString(PyExc_RuntimeError, "Internal error: invalid DHT filter index");
        return NULL;

    case -3:
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate internal variables");
        return NULL;

    case -4:
        PyErr_SetString(PyExc_RuntimeError, "Failed to converge");
        return NULL;

    case -5:
        PyErr_SetString(PyExc_ZeroDivisionError, "Internal error: division by zero");
        return NULL;

    case -6:
        PyErr_SetString(PyExc_ValueError,
                        "Internal error: wrong nzeros in function bessel_j_zero (must be >= 1)");
        return NULL;

    case -7:
        PyErr_SetString(PyExc_ValueError,
                        "Internal error: wrong n of iterations in pade sum (must be >= 1)");
        return NULL;

    case -8:
        PyErr_SetString(PyExc_ValueError, "Error: n_eval must be provided and cannot be zero");
        return NULL;

    case -9:
        PyErr_SetString(PyExc_ValueError, "Error: eps_rel must be provided and cannot be zero");
        return NULL;

    case -10:
        PyErr_SetString(PyExc_ValueError, "Error: f_max must be provided and cannot be zero");
        return NULL;

    case -11:
        PyErr_SetString(PyExc_ValueError,
                        "Error: invalid strategy name, must be one of : " LIBHANKEL_ALL_STRATEGIES
                        ".");
        return NULL;

    case -12:
        PyErr_SetString(PyExc_ValueError, "Error: x must be finite and greater than zero");
        return NULL;

    default:
        PyErr_SetString(PyExc_RuntimeError, "unknown error");
        return NULL;
    }

    PyObject *out_list = PyList_New(len_x);
    for (Py_ssize_t i = 0; i < len_x; i++) {
        PyList_SetItem(out_list, i, PyFloat_FromDouble(output[i]));
    }

    free(x);
    free(output);

    free(f_ctx->params);
    if (f_ctx->callable) {
        Py_DECREF(f_ctx->callable);
    }
    free(f_ctx);
    return out_list;
}

static PyObject *py_dht_nodes(PyObject *self, PyObject *args) {
    const char *strategy_name;
    int nu;

    if (!PyArg_ParseTuple(args, "si", &strategy_name, &nu)) {
        return NULL;
    }

    if (nu != 0 && nu != 1) {
        PyErr_SetString(PyExc_ValueError, "nu must be 0 or 1");
        return NULL;
    }

    size_t n_nodes = 0;
    const double *table = dht_node_table(strategy_name, &n_nodes);
    if (!table) {
        PyErr_SetString(PyExc_ValueError,
                        "Error: nodes are only tabulated for the fixed-abscissa filters, one of: "
                        "'DHT_Key_51', 'DHT_Key_101', 'DHT_Key_201', 'DHT_Anderson_801'.");
        return NULL;
    }

    PyObject *abscissae = PyList_New(n_nodes);
    PyObject *weights = PyList_New(n_nodes);
    if (!abscissae || !weights) {
        Py_XDECREF(abscissae);
        Py_XDECREF(weights);
        return NULL;
    }

    for (size_t i = 0; i < n_nodes; i++) {
        PyList_SetItem(abscissae, (Py_ssize_t)i, PyFloat_FromDouble(table[i * 3]));
        PyList_SetItem(weights, (Py_ssize_t)i, PyFloat_FromDouble(table[i * 3 + nu + 1]));
    }

    /* PyTuple_Pack takes its own references, so drop ours. */
    PyObject *result = PyTuple_Pack(2, abscissae, weights);
    Py_DECREF(abscissae);
    Py_DECREF(weights);
    return result;
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

// docstring for dht_nodes python api
static char dht_nodes_doc[] =
    "Return the abscissae and weights of a fixed-abscissa digital filter.\n"
    "\n"
    "The filters evaluate the transform at ``x`` as "
    "``sum(f(a[i]/x) * (a[i]/x) * w[i] / x)``, so the points at which ``f`` is "
    "needed are known before ``f`` is evaluated. Callers that already have a "
    "vectorised way to compute ``f`` can use this to gather every abscissa "
    "first, evaluate once, and combine, instead of going through "
    ":func:`hankel_transform` and its per-point callback.\n"
    "\n"
    ":param strategy_name:   The name of the filter. Only the fixed-abscissa "
    "filters are tabulated: 'DHT_Key_51', 'DHT_Key_101', 'DHT_Key_201', "
    "'DHT_Anderson_801'. The adaptive strategies choose their abscissae from "
    "values they have already seen, so they have no table.\n"
    ":type strategy_name:    str\n"
    ":param nu:              The order of the Bessel function, must be 0 or 1. "
    "Selects which column of weights is returned.\n"
    ":type nu:               int\n"
    ":returns:               ``(abscissae, weights)``, both of length equal to "
    "the filter length.\n"
    ":rtype:                 tuple[list[float], list[float]]\n";

static PyMethodDef Methods[] = {
    {"hankel_transform", py_hankel_transform, METH_VARARGS, hankel_t_doc},
    {"dht_nodes", py_dht_nodes, METH_VARARGS, dht_nodes_doc},
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef module = {PyModuleDef_HEAD_INIT, "libhankel", " ", -1, Methods};

PyMODINIT_FUNC PyInit_libhankel(void) { return PyModule_Create(&module); }
