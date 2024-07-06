#include <Python.h>

// Пример функции для основного модуля numera
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
    return PyModule_Create(&numera_module);
}