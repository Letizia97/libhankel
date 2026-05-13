#include <Python.h>
#include "form_factors.h"   
#include "libhankel.h"   

static PyObject* py_form_factor_g_dab(PyObject* self, PyObject* args) {
    double q;
    PyObject* params_obj;

    if (!PyArg_ParseTuple(args, "dO", &q, &params_obj)) {
        return NULL;
    }

    if (!PySequence_Check(params_obj)) {
        PyErr_SetString(PyExc_TypeError, "params must be a sequence");
        return NULL;
    }

    Py_ssize_t n = PySequence_Size(params_obj);

    if (n == 0) {
        PyErr_SetString(PyExc_ValueError, "params cannot be empty");
        return NULL;
    }

    if (n < 0) return NULL;  // error already set

    if (n > 50) {
        PyErr_SetString(PyExc_ValueError, "params must have at most 50 elements");
        return NULL;
    }

    // Initialize full array to 0
    double params[50] = {0};

    for (Py_ssize_t i = 0; i < n; i++) {
        PyObject* item = PySequence_GetItem(params_obj, i);
        if (!item) return NULL;

        double value = PyFloat_AsDouble(item);
        Py_DECREF(item);

        if (PyErr_Occurred()) return NULL;

        params[i] = value;
    }

    double result = form_factor_g_dab(q, params);
    return PyFloat_FromDouble(result);
}



static PyObject* py_hankel_transform(PyObject* self, PyObject* args) {

    int nu;
    PyObject *func_name, *x_obj, *params_obj, *strategy_param_obj;
    const char *strategy_name;

    // Parse Python inputs
    if (!PyArg_ParseTuple(args, "iOOOsO",
                          &nu,
                          &func_name,
                          &x_obj,
                          &params_obj,
                          &strategy_name,
                          &strategy_param_obj)) {
        return NULL;
    }

    // ---------------------------
    // Resolve function pointer
    // ---------------------------
    double (*func_ptr)(double, double (*)[50]) = NULL;

    // if (PyCallable_Check(func_name)) {
    //     // Python callback
    //     Py_XINCREF(func_name);
    //     global_py_func = func_name;

    //     func_ptr = callback_wrapper;

    // } else 
    if (PyUnicode_Check(func_name)) {
        // String → choose built-in function
        const char *name = PyUnicode_AsUTF8(func_name);

        if (strcmp(name, "g_dab") == 0) {
            func_ptr = form_factor_g_dab;
        } else if (strcmp(name, "broad_peak") == 0) {
            func_ptr = form_factor_broad_peak;
        } else if (strcmp(name, "spheres") == 0) {
            func_ptr = form_factor_sphere;
        } else {
            PyErr_SetString(PyExc_ValueError, "Unknown function name");
            return NULL;
        }

    } else {
        PyErr_SetString(PyExc_TypeError,
                        "func must be one of: g_dab, broad_peak, spheres");
        return NULL;
    }


    // ---------------------------
    // Convert x → C array
    // ---------------------------
    if (!PySequence_Check(x_obj)) {
        PyErr_SetString(PyExc_TypeError, "x must be a sequence");
        return NULL;
    }

    Py_ssize_t len_x = PySequence_Size(x_obj);
    if (len_x < 0) return NULL;

    double *x = malloc(len_x * sizeof(double));
    if (!x) {
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate x");
        return NULL;
    }

    for (Py_ssize_t i = 0; i < len_x; i++) {
        PyObject* item = PySequence_GetItem(x_obj, i);
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
    double f_params[50] = {0};

    if (!PySequence_Check(params_obj)) {
        free(x);
        PyErr_SetString(PyExc_TypeError, "params must be a sequence");
        return NULL;
    }

    Py_ssize_t n_params = PySequence_Size(params_obj);
    if (n_params > 50) {
        free(x);
        PyErr_SetString(PyExc_ValueError, "params max size = 50");
        return NULL;
    }

    for (Py_ssize_t i = 0; i < n_params; i++) {
        PyObject* item = PySequence_GetItem(params_obj, i);
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
        PyErr_SetString(PyExc_TypeError,
                        "strategy_params must be dict");
        return NULL;
    }

    strategy_params sp;
    PyObject* n_eval_obj = PyDict_GetItemString(strategy_param_obj, "n_eval");
    PyObject* eps_rel_obj = PyDict_GetItemString(strategy_param_obj, "eps_rel");
    PyObject* f_max_obj = PyDict_GetItemString(strategy_param_obj, "f_max");

    if (!n_eval_obj) {
        sp.n_eval = PyFloat_AsDouble(PyFloat_FromDouble(0.0));
    } else {
        sp.n_eval = PyFloat_AsDouble(n_eval_obj);
    }

    if (!eps_rel_obj) {
        sp.eps_rel = PyFloat_AsDouble(PyFloat_FromDouble(0.0));
    } else {
        sp.eps_rel = PyFloat_AsDouble(eps_rel_obj);
    }

    if (!f_max_obj) {
        sp.f_max = PyFloat_AsDouble(PyFloat_FromDouble(0.0));
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
    // Call your C function
    // ---------------------------
    hankel_transform(
        nu,
        func_ptr,
        x,
        &f_params,
        output,
        (size_t)len_x,
        strategy_name,
        sp
    );


    // ---------------------------
    // Convert output → Python list
    // ---------------------------
    PyObject *result = PyList_New(len_x);

    for (Py_ssize_t i = 0; i < len_x; i++) {
        PyList_SetItem(result, i,
            PyFloat_FromDouble(output[i]));
    }


    // ---------------------------
    // Cleanup
    // ---------------------------
    free(x);
    free(output);

    return result;
}

static PyMethodDef Methods[] = {
    {"form_factor_g_dab", py_form_factor_g_dab, METH_VARARGS, "Wrapped function"},
    {"hankel_transform", py_hankel_transform, METH_VARARGS, "Wrapped function"},
    {NULL, NULL, 0, NULL}
};


static struct PyModuleDef module = {
    PyModuleDef_HEAD_INIT,
    "_core",
    NULL,
    -1,
    Methods
};

PyMODINIT_FUNC PyInit__core(void) {
    return PyModule_Create(&module);
}
