#include "harmonic.h"

//прямое преобразование Фурье
/*std::vector<std::complex<double>> fft(const std::vector<double>& y) {
    size_t N = y.size();
    std::vector<std::complex<double>> Y(N);

    for (int i = 0; i < N; i++) {
        double real = 0.0;
        double imag = 0.0;
        for (int j = 0; j < N; j++) {
            real += y[j] * cos(2 * PI * i * j / N);
            imag += y[j] * sin(2 * PI * i * j / N);
        }
        Y[i] = std::complex<double>(real, imag);
    }
    return Y;
}std::vector<std::complex<double>> fft(const std::vector<double>& y) {
    size_t N = y.size();
    std::vector<std::complex<double>> Y(N);

    __m256d temp = _mm256_mul_pd(_mm256_set1_pd(2.0 / N), _mm256_set1_pd(PI));

    for (int i = 0; i < N; i++) {
        __m256d real_sum = _mm256_setzero_pd();
        __m256d imag_sum = _mm256_setzero_pd();
        __m256d i_v = _mm256_set1_pd(i);
        __m256d N_v = _mm256_set1_pd(N);

        for (int j = 0; j < N; j += 4) {
            __m256d j_v = _mm256_set_pd(j + 3, j + 2, j + 1, j);
            __m256d y_v = _mm256_set_pd(y[j + 3], y[j + 2], y[j + 1], y[j]);

            __m256d angle = _mm256_mul_pd(temp, _mm256_mul_pd(i_v, j_v));

            __m256d cos_v = _mm256_cos_pd(angle);
            __m256d sin_v = _mm256_sin_pd(angle);

            real_sum = _mm256_add_pd(real_sum, _mm256_mul_pd(y_v, cos_v));
            imag_sum = _mm256_add_pd(imag_sum, _mm256_mul_pd(y_v, sin_v));
        }

        double real[4], imag[4];
        _mm256_storeu_pd(real, real_sum);
        _mm256_storeu_pd(imag, imag_sum);
        double real_total = real[0] + real[1] + real[2] + real[3];
        double imag_total = imag[0] + imag[1] + imag[2] + imag[3];

        Y[i] = std::complex<double>(real_total, imag_total);
    }
    return Y;
}
std::vector<std::complex<double>> fft(const std::vector<double>& y) {
    int n = static_cast<int>(y.size());
    std::vector<std::complex<double>> p(n);

    if (n == 1) {
        p[0] = y[0];
        return p;
    }

    std::vector<double> a(n / 2), b(n / 2);
    for (int i = 0; i < n / 2; i++) {
        a[i] = y[2 * i];
        b[i] = y[2 * i + 1];
    }

    std::vector<std::complex<double>> A = fft(a);
    std::vector<std::complex<double>> B = fft(b);

    std::complex<double> w = 1;
    std::complex<double> wn = std::exp(-2.0 * std::complex<double>(0, 1) * PI / static_cast<double>(n));

    for (int i = 0; i < n / 2; i++) {
        p[i] = A[i] + w * B[i];
        p[i + n / 2] = A[i] - w * B[i];
        w *= wn;
    }

    return p;
}
*/


std::vector<std::complex<double>> ft(const std::vector<double>& y) {
    size_t N = y.size();
    std::vector<std::complex<double>> Y(N);

    __m256d temp = _mm256_set1_pd(2.0 * PI / N);

    for (int i = 0; i < N; i++) {
        __m256d real_sum = _mm256_setzero_pd();
        __m256d imag_sum = _mm256_setzero_pd();
        __m256d i_v = _mm256_set1_pd(static_cast<double>(i));

        for (int j = 0; j < N; j += 4) {
            __m256d j_v = _mm256_set_pd(static_cast<double>(j + 3), static_cast<double>(j + 2), static_cast<double>(j + 1), static_cast<double>(j));
            __m256d y_v = _mm256_set_pd(y[j + 3], y[j + 2], y[j + 1], y[j]);

            __m256d cos_v = _mm256_cos_pd(_mm256_mul_pd(temp, _mm256_mul_pd(i_v, j_v)));
            __m256d sin_v = _mm256_sin_pd(_mm256_mul_pd(temp, _mm256_mul_pd(i_v, j_v)));

            real_sum = _mm256_add_pd(real_sum, _mm256_mul_pd(y_v, cos_v));
            imag_sum = _mm256_add_pd(imag_sum, _mm256_mul_pd(y_v, sin_v));
        }

        alignas(32) double real[4], imag[4];
        _mm256_store_pd(real, real_sum);
        _mm256_store_pd(imag, imag_sum);

        double real_total = real[0] + real[1] + real[2] + real[3];
        double imag_total = imag[0] + imag[1] + imag[2] + imag[3];

        Y[i] = std::complex<double>(real_total, imag_total);
    }
    return Y;
}

void solve(std::complex<double>* a, int n, std::complex<double> wn) {
    if (n > 1) {
        int k = (n >> 1);
        solve(a, k, wn * wn);
        solve(a + k, k, wn * wn);
        std::complex<double> w = 1;
        for (int i = 0; i < k; i++) {
            // тут нужно быть чуть аккуратней с перезаписыванием,
            // потому что мы читаем и пишем из одного и того же массива
            std::complex<double> t = w * a[i + k];
            a[i + k] = a[i] - t;
            a[i] = a[i] + t;
            w *= wn;
        }
    }
}

std::vector<std::complex<double>> fft(const std::vector<double>& y) {
    int n = y.size();
    std::vector<std::complex<double>> a_copy(y.begin(), y.end()); // Создание копии входного массива

    const int logn = log2(n);

    for (int i = 0; i < n; i++) {
        int k = 0;
        for (int l = 0; l < logn; l++)
            k |= ((i >> l & 1) << (logn - l - 1));
        if (i < k)
            std::swap(a_copy[i], a_copy[k]);
    }

    std::complex<double> wn = std::polar(1., 1 * 2 * 3.14 / n);
    solve(&a_copy[0], n, wn); // Передача указателя на данные вектора

    return a_copy;
}
/*
std::vector<std::complex<double>> fft(const std::vector<double>& y) {
    int n = static_cast<int>(y.size());
    std::vector<std::complex<double>> p(n);

    if (n == 1) {
        p[0] = y[0];
        return p;
    }

    std::vector<double> a(n / 2), b(n / 2);
    for (int i = 0; i < n / 2; i++) {
        a[i] = y[2 * i];
        b[i] = y[2 * i + 1];
    }

    std::vector<std::complex<double>> A = fft(a);
    std::vector<std::complex<double>> B = fft(b);

    std::complex<double> w = 1;
    std::complex<double> wn = exp(-2.0 * std::complex<double>(0, 1) * PI / static_cast<double>(n));

    for (int i = 0; i < n / 2; i++) {
        p[i] = A[i] + w * B[i];
        p[i + n / 2] = A[i] - w * B[i];
        w *= wn;
    }

    return p;
}

std::vector<std::complex<double>> fft(const std::vector<double>& y) {
    int n = static_cast<int>(y.size());
    std::vector<std::complex<double>> p(n);

    if (n == 1) {
        p[0] = y[0];
        return p;
    }

    std::vector<double> a(n / 2), b(n / 2);
    for (int i = 0; i < n / 2; i++) {
        a[i] = y[2 * i];
        b[i] = y[2 * i + 1];
    }

    std::vector<std::complex<double>> A = fft(a);
    std::vector<std::complex<double>> B = fft(b);

    std::complex<double> wn = std::exp(-2.0 * std::complex<double>(0, 1) * PI / static_cast<double>(n));
    double wn_real = wn.real();
    double wn_imag = wn.imag();

    __m256d w_real = _mm256_set1_pd(1.0);
    __m256d w_imag = _mm256_set1_pd(0.0);
    __m256d wn_real_v = _mm256_set1_pd(wn_real);
    __m256d wn_imag_v = _mm256_set1_pd(0.0);

    for (int i = 0; i < n / 2; i += 4) {
        __m256d A_real = _mm256_set_pd(A[i + 3].real(), A[i + 2].real(), A[i + 1].real(), A[i].real());
        __m256d A_imag = _mm256_set_pd(A[i + 3].imag(), A[i + 2].imag(), A[i + 1].imag(), A[i].imag());
        __m256d B_real = _mm256_set_pd(B[i + 3].real(), B[i + 2].real(), B[i + 1].real(), B[i].real());
        __m256d B_imag = _mm256_set_pd(B[i + 3].imag(), B[i + 2].imag(), B[i + 1].imag(), B[i].imag());

        __m256d twiddle_real = _mm256_sub_pd(_mm256_mul_pd(w_real, B_real), _mm256_mul_pd(w_imag, B_imag));
        __m256d twiddle_imag = _mm256_add_pd(_mm256_mul_pd(w_real, B_imag), _mm256_mul_pd(w_imag, B_real));
    
        __m256d p_real = _mm256_add_pd(A_real, twiddle_real);
        __m256d p_imag = _mm256_add_pd(A_imag, twiddle_imag);
        __m256d p_half_real = _mm256_sub_pd(A_real, twiddle_real);
        __m256d p_half_imag = _mm256_sub_pd(A_imag, twiddle_imag);

        double real[4], imag[4];
        _mm256_storeu_pd(real, p_real);
        _mm256_storeu_pd(imag, p_imag);

        p[i] = std::complex<double>(real[0], imag[0]);
        p[i + 1] = std::complex<double>(real[1], imag[1]);
        p[i + 2] = std::complex<double>(real[2], imag[2]);
        p[i + 3] = std::complex<double>(real[3], imag[3]);

        _mm256_storeu_pd(real, p_half_real);
        _mm256_storeu_pd(imag, p_half_imag);

        p[i + n / 2] = std::complex<double>(real[0], imag[0]);
        p[(i + 1) + n / 2] = std::complex<double>(real[1], imag[1]);
        p[(i + 2) + n / 2] = std::complex<double>(real[2], imag[2]);
        p[(i + 3) + n / 2] = std::complex<double>(real[3], imag[3]);

        w_imag = _mm256_add_pd(_mm256_mul_pd(w_real, wn_imag_v), _mm256_mul_pd(w_imag, wn_real_v));
    }

    return p;
}*/

//обратное преобразование Фурье
std::vector<double> ift(const std::vector<std::complex<double>>& Y) {
    size_t N = Y.size();
    std::vector<double> y(N);

    double normalization = 1.0 / static_cast<double>(N);

    __m256d temp = _mm256_mul_pd(_mm256_set1_pd(2.0 * PI / N), _mm256_set1_pd(-1.0));

    for (int i = 0; i < N; ++i) {
        __m256d real_sum = _mm256_setzero_pd();
        __m256d i_v = _mm256_set1_pd(static_cast<double>(i));

        for (int j = 0; j < N; j += 4) {
            __m256d j_v = _mm256_set_pd(j + 3, j + 2, j + 1, j);
            __m256d Y_real_v = _mm256_set_pd(Y[j + 3].real(), Y[j + 2].real(), Y[j + 1].real(), Y[j].real());
            __m256d Y_imag_v = _mm256_set_pd(Y[j + 3].imag(), Y[j + 2].imag(), Y[j + 1].imag(), Y[j].imag());

            __m256d cos_v = _mm256_cos_pd(_mm256_mul_pd(temp, _mm256_mul_pd(i_v, j_v)));
            __m256d sin_v = _mm256_sin_pd(_mm256_mul_pd(temp, _mm256_mul_pd(i_v, j_v)));

            real_sum = _mm256_add_pd(real_sum, _mm256_sub_pd(
                _mm256_mul_pd(Y_real_v, cos_v), _mm256_mul_pd(Y_imag_v, sin_v)));
        }

        alignas(32) double real[4];
        _mm256_store_pd(real, real_sum);
        double real_total = real[0] + real[1] + real[2] + real[3];

        y[i] = real_total * normalization;
    }

    return y;
}

/*
std::vector<double> ifft(const std::vector<std::complex<double>>& Y) {
    size_t N = Y.size();
    std::vector<double> y(N);

    for (int i = 0; i < N; i++) {
        __m256d real = _mm256_setzero_pd();
        for (int j = 0; j < N; j++) {
            real += Y[j].real() * cos(2 * PI * i * j / N) - Y[j].imag() * sin(2 * PI * i * j / N);
        }
        y[N - i] = (real / N);
    }

    return y;
}*/

//прямое преобразование Хартли
std::vector<double> hart(const std::vector<double>& y) {
    size_t N = y.size();
    std::vector<double> H(N);

    for (int k = 0; k < N; k++) {
        double sum = 0.0;
        for (int n = 0; n < N; n++) {
            double angle = 2 * PI * n * k / N;
            sum += y[n] * (cos(angle) + sin(angle));
        }
        H[k] = sum;
    }
    return H;
}

//обратное преобразование Хартли
std::vector<double> ihart(const std::vector<double>& H) {
    size_t N = H.size();
    std::vector<double> y(N);

    for (int n = 0; n < N; n++) {
        double sum = 0.0;
        for (int k = 0; k < N; k++) {
            double angle = 2 * PI * n * k / N;
            sum += H[k] * (cos(angle) + sin(angle));
        }
        y[n] = sum / N;
    }

    return y;
}

std::vector<double> convolution(std::vector<double>& y, std::vector<double>& w) {
    size_t N1 = y.size();
    size_t N2 = w.size();
    size_t N3 = N1 + N2 - 1;

    std::vector<double> y_new(N3, 0.0);

    for (size_t i = 0; i < N3; ++i) {
        size_t j_start = (i >= N2 - 1) ? i - (N2 - 1) : 0;
        size_t j_end = (i < N1 - 1) ? i : N1 - 1;

        __m256d sum = _mm256_setzero_pd();

        for (size_t j = j_start; j <= j_end; ++j) {
            __m256d y_vec = _mm256_loadu_pd(&y[j]);
            __m256d w_vec = _mm256_set1_pd(w[i - j]);

            sum = _mm256_fmadd_pd(y_vec, w_vec, sum);
        }

        double result[4];
        _mm256_storeu_pd(result, sum);

        y_new[i] = result[0] + result[1] + result[2] + result[3];
    }

    return y_new;
}