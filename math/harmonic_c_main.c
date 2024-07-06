#include <Python.h>
#include "harmonic_c.h"

// Определение методов подмодуля
static PyMethodDef HarmonicMethods[] = {
    {"ft", ft, METH_VARARGS, "Return the input list"},
    {"fft", fft, METH_VARARGS, "Return the input list"},
    {"fftn", fftn, METH_VARARGS, "Sum the elements of a list"},
    {"hart", hart, METH_VARARGS, "Create a new bint."},
    {"ihart", ihart, METH_VARARGS, "Create a new bint."},
    {NULL, NULL, 0, NULL}
};

// Определение подмодуля
static struct PyModuleDef harmonic_module = {
    PyModuleDef_HEAD_INIT,
    "numera.harmonic",
    NULL,
    -1,
    HarmonicMethods
};

PyMODINIT_FUNC PyInit_harmonic(void) {
    return PyModule_Create(&harmonic_module);
}