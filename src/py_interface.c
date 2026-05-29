
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdlib.h>

#include "libhankel.h"
#include "form_factors.h"


double python_form_factor(
    double x,
    double *params,
    size_t n_params,
    void *user_data
) {
    PyObject *callable = (PyObject *)user_data;

    PyGILState_STATE gstate = PyGILState_Ensure();

    PyObject *args = PyTuple_New(2);
    PyTuple_SetItem(args, 0, PyFloat_FromDouble(x));

    PyObject *list = PyList_New(n_params);
    for (size_t i = 0; i < n_params; i++) {
        PyList_SetItem(list, i, PyFloat_FromDouble(params[i]));
    }

    PyTuple_SetItem(args, 1, list);

    PyObject *result = PyObject_CallObject(callable, args);
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

static PyObject* py_hankel_transform(PyObject *self, PyObject *args) {
    int nu;
    PyObject *f_obj, *x_obj, *params_obj, *strategy_param_obj;
    const char *strategy_name;

    if (!PyArg_ParseTuple(
            args,
            "iOOOsO",
            &nu,
            &f_obj,
            &x_obj,
            &params_obj,
            &strategy_name,
            &strategy_param_obj)) 
        return NULL;


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
    if (!PySequence_Check(params_obj)) {
        free(x);
        PyErr_SetString(PyExc_TypeError, "params must be a sequence");
        return NULL;
    }

    Py_ssize_t n_params = PySequence_Size(params_obj);
    double *f_params = malloc(n_params * sizeof(double));

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
        PyErr_SetString(
            PyExc_TypeError,
            "strategy_params must be dict"
        );
        return NULL;
    }

    strategy_params sp;
    PyObject* n_eval_obj = PyDict_GetItemString(strategy_param_obj, "n_eval");
    PyObject* eps_rel_obj = PyDict_GetItemString(strategy_param_obj, "eps_rel");
    PyObject* f_max_obj = PyDict_GetItemString(strategy_param_obj, "f_max");

    if (!n_eval_obj) {
        sp.n_eval = 0.0;
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
    void *user_data = NULL;

    /* check for Python callable */
    if (PyCallable_Check(f_obj)) {
        f_ptr = python_form_factor;
        user_data = f_obj;
        Py_INCREF(f_obj);

    /* or whether the user wants to use built-in form factor */
    } else if (PyUnicode_Check(f_obj)) {
        const char *name = PyUnicode_AsUTF8(f_obj);
        f_ptr = get_form_factor_by_name(name);

        if (!f_ptr) {
            free(x); free(f_params); free(output);
            PyErr_SetString(PyExc_ValueError, "Unknown function name");
            return NULL;
        }

    } else {
        free(x); free(f_params); free(output);
        PyErr_SetString(PyExc_TypeError, "f must be callable or string");
        return NULL;
    }

    hankel_transform(
        nu,
        f_ptr,
        x,
        len_x,
        f_params,
        n_params,
        output,
        strategy_name,
        sp,
        user_data
    );

    PyObject *out_list = PyList_New(len_x);
    for (Py_ssize_t i = 0; i < len_x; i++) {
        PyList_SetItem(out_list, i, PyFloat_FromDouble(output[i]));
    }

    free(x);
    free(f_params);
    free(output);

    if (PyCallable_Check(f_obj)) {
        Py_DECREF(f_obj);
    }

    return out_list;
}

static PyMethodDef Methods[] = {
    {"hankel_transform", py_hankel_transform, METH_VARARGS, "Hankel transform"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef module = {
    PyModuleDef_HEAD_INIT,
    "libhankel",
    NULL,
    -1,
    Methods
};

PyMODINIT_FUNC PyInit_libhankel(void) {
    return PyModule_Create(&module);
}
