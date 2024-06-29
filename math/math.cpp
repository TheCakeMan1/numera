#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "interpol.h"
#include "parser.h"
#include <pybind11/complex.h>
#include <cmath>
#include <stack>

namespace py = pybind11;
double PI = 3.14159265358979323846264338327950288419716939937510;

/*py::object find_index(const std::vector<py::object>& items, const py::object& value) {
    for (int i = 0; i < items.size(); i++) {
        if (items[i] == value) {
            return py::cast(i);
        }
    }   
    return py::none();
}*/


std::vector<double> vec_swap(const std::vector<double>& y) {
    std::vector<double> temp;
    for (auto i : y) {
        temp.push_back(i);
    }
    return temp;
}

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
        double imag = 0.0;
        for (int j = 0; j < N; j++) {
            real += Y[j].real() * cos(2 * PI * i * j / N) - Y[j].imag() * sin(2 * PI * i * j / N);
            imag += Y[j].real() * sin(2 * PI * i * j / N) + Y[j].imag() * cos(2 * PI * i * j / N);
        }
        y[i] = (real / N) + (imag / N);
    }
    return vec_swap(y);
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

int add(int a, int b) {
    return a + b;
}

int subtract(int a, int b) {
    return a - b;
}

PYBIND11_MODULE(numera, m) {
    m.doc() = ""; // optional module docstring

    m.attr("pi") = py::float_(PI);

    // Define submodules
    py::module submodule = m.def_submodule("harmonic", "");
    submodule.def("fft", &fft, "");
    submodule.def("ifft", &ifft, "");
    submodule.def("hart", &hart, "");
    submodule.def("ihart", &ihart, "");

    // Add functions to submodule
    //submodule.def("add", &add, "Addition function");
    //submodule.def("subtract", &subtract, "Subtraction function");
	//submodule.def("parser", &parser, "Subtraction function");
    //submodule.def("find_index", &find_index, "find_index function");
}