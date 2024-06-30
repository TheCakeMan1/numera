#include "harmonic.h"

//прямое преобразование Фурье
std::vector<std::complex<double>> fft(const std::vector<double>& y) {
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
}

//обратное преобразование Фурье
std::vector<double> ifft(const std::vector<std::complex<double>>& Y) {
    size_t N = Y.size();
    std::vector<double> y(N);

    for (int i = 0; i < N; i++) {
        double real = 0.0;
        for (int j = 0; j < N; j++) {
            real += Y[j].real() * cos(2 * PI * i * j / N) - Y[j].imag() * sin(2 * PI * i * j / N);
        }
        y[N - i] = (real / N);
    }

    return y;
}

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