#ifndef NUI_H
#define NUI_H
#include <immintrin.h>
#include <corecrt_math.h>
#include <ostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

/*class bint {
private:
    std::vector<__m256i> digits;

public:
    // Конструкторы
    bint() {
        digits.push_back(_mm256_setzero_si256());
    }

    bint(int a) {
        std::string value_str = std::to_string(a);
        size_t len = value_str.size();
        size_t num_avx_elements = (len + 7) / 8;  // Each __m256i can hold 8 int values (256 bits / 32 bits per int)

        digits.resize(num_avx_elements, _mm256_setzero_si256());

        for (size_t i = 0; i < len; ++i) {
            if (value_str[len - 1 - i] >= '0' && value_str[len - 1 - i] <= '9') {
                reinterpret_cast<int*>(&digits[i / 8])[i % 8] = value_str[len - 1 - i] - '0';
            }
        }

        // Remove leading zero AVX elements
        while (digits.size() > 1 && _mm256_testz_si256(digits.back(), digits.back()) == 1) {
            digits.pop_back();
        }
    }

    // Оператор вывода
    friend std::ostream& operator<<(std::ostream& os, const bint& bigint) {
        bool leading_zero = true;
        for (auto it = bigint.digits.rbegin(); it != bigint.digits.rend(); ++it) {
            for (int i = 7; i >= 0; --i) {
                int digit = reinterpret_cast<const int*>(it->m256i_u32)[i];
                if (digit != 0) {
                    leading_zero = false;
                }
                if (!leading_zero || i == 0) {
                    os << digit;
                }
            }
        }
        return os;
    }

    // Операция сложения
    bint operator+(const bint& other) const {
        bint result;
        size_t max_len = std::max(digits.size(), other.digits.size());
        result.digits.resize(max_len + 1, _mm256_setzero_si256());

        __m256i carry = _mm256_setzero_si256();

        for (size_t i = 0; i < max_len; ++i) {
            __m256i a = (i < digits.size()) ? digits[i] : _mm256_setzero_si256();
            __m256i b = (i < other.digits.size()) ? other.digits[i] : _mm256_setzero_si256();

            __m256i sum = _mm256_add_epi32(a, b);
            sum = _mm256_add_epi32(sum, carry);

            __m256i carry_mask = _mm256_cmpgt_epi32(sum, _mm256_set1_epi32(9));
            carry = _mm256_and_si256(carry_mask, _mm256_set1_epi32(1));

            sum = _mm256_sub_epi32(sum, _mm256_and_si256(carry_mask, _mm256_set1_epi32(10)));

            result.digits[i] = sum;
        }

        if (_mm256_testz_si256(carry, carry) == 0) {
            result.digits[max_len] = carry;
        }
        else {
            result.digits.pop_back();
        }

        return result;
    }

    // Операция умножения
    bint operator*(const bint& other) const {
        size_t new_size = digits.size() + other.digits.size();
        std::vector<__m256i> result_digits(new_size, _mm256_setzero_si256());

        for (size_t i = 0; i < digits.size(); ++i) {
            __m256i carry = _mm256_setzero_si256();
            for (size_t j = 0; j < other.digits.size(); ++j) {
                __m256i prod = _mm256_setzero_si256();
                for (int k = 0; k < 8; ++k) {
                    int* result_ptr = reinterpret_cast<int*>(&result_digits[i + j]);
                    const int* digit_ptr = reinterpret_cast<const int*>(&digits[i]);
                    const int* other_digit_ptr = reinterpret_cast<const int*>(&other.digits[j]);
                    int* carry_ptr = reinterpret_cast<int*>(&carry);

                    long long current = result_ptr[k] + digit_ptr[k] * 1LL * other_digit_ptr[k] + carry_ptr[k];
                    reinterpret_cast<int*>(&prod)[k] = current % 10;
                    carry_ptr[k] = current / 10;
                }
                result_digits[i + j] = prod;
            }
            // Propagate carry
            for (size_t k = other.digits.size(); k < result_digits.size(); ++k) {
                if (_mm256_testz_si256(carry, carry) == 1) break;
                __m256i current = _mm256_add_epi32(result_digits[i + k], carry);
                __m256i carry_mask = _mm256_cmpgt_epi32(current, _mm256_set1_epi32(9));
                carry = _mm256_and_si256(carry_mask, _mm256_set1_epi32(1));
                current = _mm256_sub_epi32(current, _mm256_and_si256(carry_mask, _mm256_set1_epi32(10)));
                result_digits[i + k] = current;
            }
        }

        while (result_digits.size() > 1 && _mm256_testz_si256(result_digits.back(), result_digits.back()) == 1) {
            result_digits.pop_back();
        }

        bint result;
        result.digits = result_digits;
        return result;
    }
};*/

class nui {
protected:
	int a;
public:
	nui() = default;
	nui(int b) {
		this->a = static_cast<int>(b);
	}
	//ni(double b) {
	//	this->a = round(b);
	//}
	nui(const nui& b) {
		this->a = b.a;
	}
	int toInt() {
		return a;
	}
	float toFloat() {
		return a;
	}
	//nf toNf();
	nui operator+(const nui& b) const {
		int result;
		__m256i a_v = _mm256_set1_epi64x(a);
		__m256i b_v = _mm256_set1_epi64x(b.a);
		__m256i c_v = _mm256_add_epi64(a_v, b_v);
		_mm256_storeu_si256((__m256i*) & result, c_v);
		return nui(result);
	}
	nui operator+(int b) const {
		double result[1];
		__m256d a_v = _mm256_set1_pd(a);
		__m256d b_v = _mm256_set1_pd(b);
		__m256d c_v = _mm256_add_pd(a_v, b_v);
		_mm256_storeu_pd(result, c_v);
		return nui(result[0]);
	}
	nui operator+(double b) const {
		double result[1];
		__m256d a_v = _mm256_set1_pd(a);
		__m256d b_v = _mm256_set1_pd(b);
		__m256d c_v = _mm256_add_pd(a_v, b_v);
		_mm256_storeu_pd(result, c_v);
		return nui(result[0]);
	}
	//ni operator+(const nf& b) const;
	nui operator-(const nui& b) const {
		double result[1];
		__m256d a_v = _mm256_set1_pd(a);
		__m256d b_v = _mm256_set1_pd(b.a);
		__m256d c_v = _mm256_sub_pd(a_v, b_v);
		_mm256_storeu_pd(result, c_v);
		return nui(result[0]);
	}

	nui operator*(const nui& b) const {
		double result[1];
		__m256d a_v = _mm256_set1_pd(a);
		__m256d b_v = _mm256_set1_pd(b.a);
		__m256d c_v = _mm256_mul_pd(a_v, b_v);
		_mm256_storeu_pd(result, c_v);
		return nui(result[0]);
	}

	std::string print() const {
		std::ostringstream oss;
		oss << a;
		return oss.str();
	}

	friend std::ostream& operator<<(std::ostream& os, const nui& val);
	//friend std::ostream& operator<<(std::ostream& os, const ni per);
};

#endif