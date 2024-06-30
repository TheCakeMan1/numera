#ifndef FUNC_H
#define FUNC_H
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include "val.h"
#include <immintrin.h>
#include "vector.h"

array gaus_func(int& n, double& a, double& b, double& c, double& x_start, double& x_end);

#endif
