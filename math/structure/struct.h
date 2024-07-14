#ifndef STRYCT_C_H
#define STRYCT_C_H

#include <Python.h>
#include <math.h>

typedef struct {
    long double real;
    long double imag;
} Complex;

Complex create_complex(long double real, long double imag);
Complex mul_comlex(Complex x, Complex y);
Complex mul_comint(Complex x, long double y);
Complex sub_comlex(Complex x, Complex y);
Complex sub_comint(Complex x, long double y);
Complex add_comlex(Complex x, Complex y);
Complex add_comint(Complex x, long double y);
Complex polar(long double magnitude, long double angle);
Complex pycomlex_ascomlex(PyObject* x);

#endif
