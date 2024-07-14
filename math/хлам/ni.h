#ifndef NI_H
#define NI_H
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <immintrin.h>
#include <ostream>
#include "nf.h"

class ni {
private:
	int64_t* a;
public:
	ni() = default;
	ni(long long b) {
		this->a = new int64_t(static_cast<int64_t>(b));
	}
	//ni(double b) {
	//	this->a = round(b);
	//}
	ni(const ni& b) {
		this->a = b.a;
	}
	int64_t toInt() {
		return *a;
	}
	double toFloat() {
		return *a;
	}
	//nf toNf();
	ni operator+(const ni& b) const {
		int result;
		__m256i a_v = _mm256_set1_epi64x(*a);
		__m256i b_v = _mm256_set1_epi64x(*b.a);
		__m256i c_v = _mm256_add_epi64(a_v, b_v);
		_mm256_storeu_si256((__m256i*) & result, c_v);
		return ni(result);
	}
	ni operator+(int b) const {
		double result[1];
		__m256d a_v = _mm256_set1_pd(*a);
		__m256d b_v = _mm256_set1_pd(b);
		__m256d c_v = _mm256_add_pd(a_v, b_v);
		_mm256_storeu_pd(result, c_v);
		return ni(result[0]);
	}
	ni operator+(double b) const {
		double result[1];
		__m256d a_v = _mm256_set1_pd(*a);
		__m256d b_v = _mm256_set1_pd(b);
		__m256d c_v = _mm256_add_pd(a_v, b_v);
		_mm256_storeu_pd(result, c_v);
		return ni(result[0]);
	}
	//ni operator+(const nf& b) const;
	ni operator-(const ni& b) const {
		double result[1];
		__m256d a_v = _mm256_set1_pd(*a);
		__m256d b_v = _mm256_set1_pd(*b.a);
		__m256d c_v = _mm256_sub_pd(a_v, b_v);
		_mm256_storeu_pd(result, c_v);
		return ni(result[0]);
	}

	ni operator*(const ni& b) const {
		double result[1];
		__m256d a_v = _mm256_set1_pd(*a);
		__m256d b_v = _mm256_set1_pd(*b.a);
		__m256d c_v = _mm256_mul_pd(a_v, b_v);
		_mm256_storeu_pd(result, c_v);
		return ni(result[0]);
	}

	std::string print() const {
		std::ostringstream oss;
		oss << *a;
		return oss.str();
	}

	friend std::ostream& operator<<(std::ostream& os, const ni& val);
	//friend std::ostream& operator<<(std::ostream& os, const ni per);
};

#endif