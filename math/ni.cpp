#include "ni.h"

std::ostream& operator<<(std::ostream& os, const ni& per) {
    os << per.a;
    return os;
}

//nf ni::toNf() {
//    return nf(a);
//}
/*
ni ni::operator+(const ni& b) const {
    double result[1];
    __m256d a_v = _mm256_set1_pd(a);
    __m256d b_v = _mm256_set1_pd(b.a);
    __m256d c_v = _mm256_add_pd(a_v, b_v);
    _mm256_storeu_pd(result, c_v);
    return ni(result[0]);
}

ni ni::operator+(int b) const {
    double result[1];
    __m256d a_v = _mm256_set1_pd(a);
    __m256d b_v = _mm256_set1_pd(b);
    __m256d c_v = _mm256_add_pd(a_v, b_v);
    _mm256_storeu_pd(result, c_v);
    return ni(result[0]);
}

ni ni::operator+(double b) const {
    double result[1];
    __m256d a_v = _mm256_set1_pd(a);
    __m256d b_v = _mm256_set1_pd(b);
    __m256d c_v = _mm256_add_pd(a_v, b_v);
    _mm256_storeu_pd(result, c_v);
    return ni(result[0]);
}

ni ni::operator+(const nf& b) const {
    double result[1];
    __m256d a_v = _mm256_set1_pd(a);
    __m256d b_v = _mm256_set1_pd(b.toFloat());
    __m256d c_v = _mm256_add_pd(a_v, b_v);
    _mm256_storeu_pd(result, c_v);
    return ni(result[0]);
}

ni ni::operator-(const ni& b) const {
    double result[1];
    __m256d a_v = _mm256_set1_pd(a);
    __m256d b_v = _mm256_set1_pd(b.a);
    __m256d c_v = _mm256_sub_pd(a_v, b_v);
    _mm256_storeu_pd(result, c_v);
    return ni(result[0]);
}

ni ni::operator*(const ni& b) const {
    double result[1];
    __m256d a_v = _mm256_set1_pd(a);
    __m256d b_v = _mm256_set1_pd(b.a);
    __m256d c_v = _mm256_mul_pd(a_v, b_v);
    _mm256_storeu_pd(result, c_v);
    return ni(result[0]);
}*/

//std::ostream& operator<<(std::ostream& os, const ni per) {
 //   os << per.a;
 //   return os;
//}