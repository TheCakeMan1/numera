#ifndef HARMONIC_H
#define HARMONIC_H
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include "val.h"
#include "vector.h"

std::vector<std::complex<double>> fft(const std::vector<double>& y);
std::vector<double> ifft(const std::vector<std::complex<double>>& Y);
std::vector<double> hart(const std::vector<double>& y);
std::vector<double> ihart(const std::vector<double>& H);

#endif // HARMONIC_H