#include <Python.h>
#include "matan.h"

static PyObject* numera_hello(PyObject* self, PyObject* args) {
    return Py_BuildValue("s", "Hello from numera");
}

static PyMethodDef NumeraMethods[] = {
    {"hello", numera_hello, METH_VARARGS, "Return a greeting."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef numera_module = {
    PyModuleDef_HEAD_INIT,
    "numera",
    NULL,
    -1,
    NumeraMethods
};

PyMODINIT_FUNC PyInit_numera(void) {
    PyObject* module = PyModule_Create(&numera_module);
    if (!module) return NULL;

    // Инициализация подмодуля
    PyObject* harmonic_module = PyInit_numera_harmonic();
    if (!harmonic_module) {
        Py_DECREF(module);
        return NULL;
    }

    // Добавление подмодуля в основной модуль
    PyModule_AddObject(module, "harmonic", harmonic_module);

    return module;
}