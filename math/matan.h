#ifndef MATAN_C_H
#define MATAN_C_H

#include "harmonic_c_main.h"

static PyMethodDef HarmonicMethods[] = {
    {"ft", ft, METH_VARARGS, "Return a greeting."},
    {"fft", fft, METH_VARARGS, "Return a greeting."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef harmonic_module = {
    PyModuleDef_HEAD_INIT,
    "numera.harmonic",
    NULL,
    -1,
    HarmonicMethods
};

PyMODINIT_FUNC PyInit_numera_harmonic(void) {
    return PyModule_Create(&harmonic_module);
}

#endif