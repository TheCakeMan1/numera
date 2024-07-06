#include "harmonic_c.h"
#include <math.h>

#define PI 3.14159265358979323846

// Определение функции ft
PyObject* ft(PyObject* self, PyObject* args) {
    PyObject* list;

    if (!PyArg_ParseTuple(args, "O", &list))
        return NULL;

    Py_ssize_t N = PyList_Size(list);
    PyObject* output_list = PyList_New(N);

    for (Py_ssize_t k = 0; k < N; k++) {
        PyObject* real = PyFloat_FromDouble(0.0);
        PyObject* imag = PyFloat_FromDouble(0.0);

        if (real == NULL || imag == NULL) {
            Py_XDECREF(real);
            Py_XDECREF(imag);
            Py_DECREF(output_list);
            return NULL;
        }

        for (Py_ssize_t n = 0; n < N; n++) {
            PyObject* temp = PyList_GetItem(list, n);
            double angle = 2.0 * PI * k * n / N;
            double cosine = cos(angle);
            double sine = sin(angle);

            PyObject* cos_val = PyFloat_FromDouble(cosine);
            PyObject* sin_val = PyFloat_FromDouble(sine);

            if (cos_val == NULL || sin_val == NULL) {
                Py_DECREF(real);
                Py_DECREF(imag);
                Py_XDECREF(cos_val);
                Py_XDECREF(sin_val);
                Py_DECREF(output_list);
                return NULL;
            }

            PyObject* temp_cos = PyNumber_Multiply(temp, cos_val);
            PyObject* temp_sin = PyNumber_Multiply(temp, sin_val);

            Py_DECREF(cos_val);
            Py_DECREF(sin_val);

            if (temp_cos == NULL || temp_sin == NULL) {
                Py_DECREF(real);
                Py_DECREF(imag);
                Py_XDECREF(temp_cos);
                Py_XDECREF(temp_sin);
                Py_DECREF(output_list);
                return NULL;
            }

            PyObject* new_real = PyNumber_Add(real, temp_cos);
            PyObject* new_imag = PyNumber_Add(imag, temp_sin);

            Py_DECREF(real);
            Py_DECREF(imag);
            Py_DECREF(temp_cos);
            Py_DECREF(temp_sin);

            if (new_real == NULL || new_imag == NULL) {
                Py_XDECREF(new_real);
                Py_XDECREF(new_imag);
                Py_DECREF(output_list);
                return NULL;
            }

            real = new_real;
            imag = new_imag;
        }

        PyObject* complex_num = PyComplex_FromDoubles(
            PyFloat_AsDouble(real),
            PyFloat_AsDouble(imag)
        );

        Py_DECREF(real);
        Py_DECREF(imag);

        if (complex_num == NULL) {
            Py_DECREF(output_list);
            return NULL;
        }

        PyList_SetItem(output_list, k, complex_num);
    }

    return output_list;
}
/*
void solve(PyObject* a, Py_ssize_t n, PyObject* wn) {
    if (n > 1) {
        Py_ssize_t k = n >> 1;
        solve(a, k, PyNumber_Multiply(wn, wn));
        solve((PyObject*)((PyListObject*)a)->ob_item + k, k, PyNumber_Multiply(wn, wn));
        PyObject* w = PyComplex_FromDoubles(1.0, 0.0);
        for (int i = 0; i < k; i++) {
            PyObject* t = PyNumber_Multiply(w, PyList_GetItem(a, i + k));
            PyList_SetItem(a, i + k, PyNumber_Subtract(PyList_GetItem(a, i), t));
            PyList_SetItem(a, i, PyNumber_Add(PyList_GetItem(a, i), t));
            Py_DECREF(t);
            w = PyNumber_Multiply(w, wn);
        }
        Py_DECREF(w);
    }
}

PyObject* copy_list(PyObject* original_list) {
    Py_ssize_t size = PyList_Size(original_list);
    if (size == -1) {
        // Обработка ошибки, если original_list не является списком
        return NULL;
    }

    PyObject* copied_list = PyList_New(size);
    if (copied_list == NULL) {
        // Обработка ошибки, если не удалось создать новый список
        return NULL;
    }

    for (Py_ssize_t i = 0; i < size; ++i) {
        PyObject* item = PyList_GetItem(original_list, i);
        Py_INCREF(item);  // Увеличиваем счетчик ссылок для каждого элемента
        if (PyList_SetItem(copied_list, i, item) == -1) {
            // Обработка ошибки при установке элемента в скопированный список
            Py_DECREF(copied_list);
            return NULL;
        }
    }

    return copied_list;
}

PyObject* fft(PyObject* self, PyObject* args) {
    PyObject* list_obj;

    if (!PyArg_ParseTuple(args, "O", &list_obj))
        return NULL;
    
    PyObject* copied_list = copy_list(list_obj);
    Py_ssize_t n = PyList_Size(list_obj);
    Py_ssize_t logn = (Py_ssize_t)log2(n);

    for (Py_ssize_t i = 0; i < n; i++) {
        Py_ssize_t k = 0;
        for (Py_ssize_t l = 0; l < logn; l++)
            k |= ((i >> l) & 1) << (logn - l - 1);
        if (i < k) {
            PyObject* temp = PyList_GetItem(list_obj, i);
            PyList_SetItem(copied_list, i, PyList_GetItem(copied_list, k));
            PyList_SetItem(copied_list, k, temp);
        }
    }

    PyObject* wn = PyComplex_FromDoubles(cos(1 * 2.0 * PI / n), sin(1 * 2.0 * PI / n));
    solve(copied_list, n, wn);
    Py_DECREF(wn);
    return copied_list;
}*/

PyObject* fftn(PyObject* self, PyObject* args) {
    PyObject* list;

    if (!PyArg_ParseTuple(args, "O", &list))
        return NULL;
    
    Py_ssize_t N = PyList_Size(list);
    if (N <= 1) {
        Py_INCREF(list);
        return list;
    }

    PyObject* even = PyList_New(N / 2);
    PyObject* odd = PyList_New(N / 2);

    for (Py_ssize_t i = 0; i < N / 2; i++) {
        PyList_SetItem(even, i, PyList_GetItem(list, i * 2));
        PyList_SetItem(odd, i, PyList_GetItem(list, i * 2 + 1));
    }

    PyObject* even_fft_args = Py_BuildValue("(O)", even);
    PyObject* odd_fft_args = Py_BuildValue("(O)", odd);

    PyObject* even_fft = fftn(self, even_fft_args);
    PyObject* odd_fft = fftn(self, odd_fft_args);

    Py_DECREF(even_fft_args);
    Py_DECREF(odd_fft_args);

    PyObject* output_list = PyList_New(N);

    for (Py_ssize_t k = 0; k < N / 2; k++) {
        PyComplexObject* t;
        t = PyComplex_FromDoubles(cos(-2.0 * PI * k / N) * PyComplex_RealAsDouble(PyList_GetItem(odd_fft, k)) - sin(-2.0 * PI * k / N) * PyComplex_ImagAsDouble(PyList_GetItem(odd_fft, k)), sin(-2.0 * PI * k / N) * PyComplex_RealAsDouble(PyList_GetItem(odd_fft, k)) + cos(-2.0 * PI * k / N) * PyComplex_ImagAsDouble(PyList_GetItem(odd_fft, k)));

        PyObject* even_k = PyList_GetItem(even_fft, k);
        PyObject* even_k_plus_t = PyNumber_Add(even_k, t);
        PyObject* even_k_minus_t = PyNumber_Subtract(even_k, t);

        PyList_SetItem(output_list, k, even_k_plus_t);
        PyList_SetItem(output_list, k + N / 2, even_k_minus_t);

        Py_DECREF(t);
        Py_DECREF(even_k_plus_t);
        Py_DECREF(even_k_minus_t);
    }

    Py_DECREF(even_fft);
    Py_DECREF(odd_fft);
    Py_DECREF(even);
    Py_DECREF(odd);

    return output_list;
}

typedef struct {
    double real;
    double imag;
} Complex;

void fft(Complex* X, int n) {
    if (n <= 1) return;

    Complex* even = (Complex*)malloc(n / 2 * sizeof(Complex));
    Complex* odd = (Complex*)malloc(n / 2 * sizeof(Complex));

    for (int i = 0; i < n / 2; i++) {
        even[i] = X[i * 2];
        odd[i] = X[i * 2 + 1];
    }

    fft(even, n / 2);
    fft(odd, n / 2);

    for (int k = 0; k < n / 2; k++) {
        double t = -2 * PI * k / n;
        Complex e = { cos(t), sin(t) };
        Complex temp = { e.real * odd[k].real - e.imag * odd[k].imag,
                        e.real * odd[k].imag + e.imag * odd[k].real };

        X[k].real = even[k].real + temp.real;
        X[k].imag = even[k].imag + temp.imag;
        X[k + n / 2].real = even[k].real - temp.real;
        X[k + n / 2].imag = even[k].imag - temp.imag;
    }

    free(even);
    free(odd);
}

static PyObject* py_fft(PyObject* self, PyObject* args) {
    PyObject* input_list;
    if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &input_list)) {
        return NULL;
    }

    Py_ssize_t n = PyList_Size(input_list);
    if (n & (n - 1)) {
        PyErr_SetString(PyExc_ValueError, "Size of input list must be a power of 2");
        return NULL;
    }

    Complex* data = (Complex*)malloc(n * sizeof(Complex));
    for (Py_ssize_t i = 0; i < n; i++) {
        PyObject* item = PyList_GetItem(input_list, i);
        if (PyComplex_Check(item)) {
            data[i].real = PyComplex_RealAsDouble(item);
            data[i].imag = PyComplex_ImagAsDouble(item);
        }
        else if (PyFloat_Check(item)) {
            data[i].real = PyFloat_AsDouble(item);
            data[i].imag = 0.0;
        }
        else if (PyLong_Check(item)) {
            data[i].real = PyLong_AsDouble(item);
            data[i].imag = 0.0;
        }
        else {
            free(data);
            PyErr_SetString(PyExc_TypeError, "List elements must be numbers");
            return NULL;
        }
    }

    fft(data, n);

    PyObject* result_list = PyList_New(n);
    for (Py_ssize_t i = 0; i < n; i++) {
        PyObject* complex_num = PyComplex_FromDoubles(data[i].real, data[i].imag);
        PyList_SetItem(result_list, i, complex_num);
    }

    free(data);
    return result_list;
}