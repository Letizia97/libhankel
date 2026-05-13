#include <Python.h>
#include "form_factors.h"   
 
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

static PyMethodDef Methods[] = {
    {"form_factor_g_dab", py_form_factor_g_dab, METH_VARARGS, "Wrapped function"},
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