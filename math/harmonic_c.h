#ifndef HARMONIC_C_H
#define HARMONIC_C_H

#include <Python.h>
#include <complex.h>
#include <math.h>

// Объявление функции ft
PyObject* ft(PyObject* self, PyObject* args);
PyObject* py_fft(PyObject* self, PyObject* args);
PyObject* fftn(PyObject* self, PyObject* args);

#endif // HARMONIC_C_H

/*
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include "val.h"
#include <immintrin.h>
#include "vector.h"

std::vector<std::complex<double>> fft(const std::vector<double>& y);
std::vector<double> ift(const std::vector<std::complex<double>>& Y);
std::vector<double> hart(const std::vector<double>& y);
std::vector<double> ihart(const std::vector<double>& H);
std::vector<double> convolution(std::vector<double>& y, std::vector<double>& w);
void solve(std::complex<double>* a, int n, std::complex<double> wn);
*/
