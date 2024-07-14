#include "struct.h"

Complex create_complex(long double real, long double imag) {
	Complex temp;
	temp.real = real;
	temp.imag = imag;
	return temp;
}

Complex mul_comlex(Complex x, Complex y) {
	Complex temp;
	temp.real = x.real * y.real;
	temp.imag = x.imag * y.imag;
	return temp;
}

Complex mul_comint(Complex x, long double y) {
	Complex temp;
	temp.real = x.real * y;
	temp.imag = x.imag * y;
	return temp;
}

Complex sub_comint(Complex x, long double y) {
	Complex temp;
	temp.real = x.real - y;
	temp.imag = x.imag;
	return temp;
}

Complex sub_comlex(Complex x, Complex y) {
	Complex temp;
	temp.real = x.real - y.real;
	temp.imag = x.imag - y.imag;
	return temp;
}

Complex add_comint(Complex x, long double y) {
	Complex temp;
	temp.real = x.real + y;
	temp.imag = x.imag;
	return temp;
}

Complex add_comlex(Complex x, Complex y) {
	Complex temp;
	temp.real = x.real + y.real;
	temp.imag = x.imag + y.imag;
	return temp;
}

Complex polar(long double magnitude, long double angle) {
	Complex result;
	result.real = magnitude * cos(angle);
	result.imag = magnitude * sin(angle);
	return result;
}

Complex pycomlex_ascomlex(PyObject* x) {
	Complex temp;
	temp.real = PyComplex_RealAsDouble(x);
	temp.imag = PyComplex_ImagAsDouble(x);
	return temp;
}