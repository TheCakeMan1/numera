#ifndef PER_H
#define PER_H
#include <immintrin.h>
#include <corecrt_math.h>
#include <ostream>
#include "ni.h"

class nf {
private:
	double a;
public:
	nf() = default;
	nf(double b) {
		this->a = b;
	}
	nf(int b) {
		this->a = b;
	}
	nf(const nf& b) {
		this->a = b.a;
	}

	double toFloat() const {
		return a;
	}
};
#endif