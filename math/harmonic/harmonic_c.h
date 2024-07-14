#ifndef HARMONIC_C_H
#define HARMONIC_C_H

#include <Python.h>
#include <math.h>
#include "../structure/struct.h"
#include <immintrin.h>

PyObject* ft(PyObject* self, PyObject* args);
PyObject* fft(PyObject* self, PyObject* args);
PyObject* fftn(PyObject* self, PyObject* args);
PyObject* hart(PyObject* self, PyObject* args);
PyObject* ihart(PyObject* self, PyObject* args);

#endif // HARMONIC_C_H
