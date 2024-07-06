#include "harmonic_c.h"
#include <math.h>

#define PI 3.14159265358979323846
#define DPI 6.28318530717958647692

//#define __AVX_BUILD__ = 0

#ifdef __AVX_BUILD__

void solve(Complex* a, int n, Complex wn) {
    if (n > 1) {
        int k = (n >> 1);
        solve(a, k, mul_comlex(wn, wn));
        solve(a + k, k, mul_comlex(wn, wn));
        Complex w = create_complex(1.0, 0.0);

        for (int i = 0; i < k; i += 2) {
            __m256d w_re = _mm256_set1_pd(w.real);
            __m256d w_im = _mm256_set1_pd(w.imag);

            __m256d a_re = _mm256_set_pd(a[i + k + 1].real, a[i + k].real, a[i + 1].real, a[i].real);
            __m256d a_im = _mm256_set_pd(a[i + k + 1].imag, a[i + k].imag, a[i + 1].imag, a[i].imag);

            __m256d t_re = _mm256_sub_pd(_mm256_mul_pd(w_re, a_re), _mm256_mul_pd(w_im, a_im));
            __m256d t_im = _mm256_add_pd(_mm256_mul_pd(w_re, a_im), _mm256_mul_pd(w_im, a_re));

            _mm256_storeu_pd((double*)(&a[i + k]), _mm256_sub_pd(a_re, t_re));
            _mm256_storeu_pd((double*)(&a[i]), _mm256_add_pd(a_re, t_re));

            w = mul_comlex(w, wn);
        }
    }
}

#else

PyObject* ft(PyObject* self, PyObject* args) {
    //создаём исходный список
    PyObject* list;

    //проверка списка на существование
    if (!PyArg_ParseTuple(args, "O", &list))
        return NULL;

    Py_ssize_t N = PyList_Size(list); //размер входного и выходного массива
    PyObject* output_list = PyList_New(N); //массив преобразования фурье

    //основной цикл преобразования Фурье
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

void solve(Complex* a, int n, Complex wn) {
    if (n > 1) {
        int k = (n >> 1);
        solve(a, k, mul_comlex(wn, wn));
        solve(a + k, k, mul_comlex(wn, wn));
        Complex w;
        w.real = 1.0;
        w.imag = 0.0;
        for (Py_ssize_t i = 0; i < k; i++) {
            Complex t = mul_comlex(w, a[i + k]);
            a[i + k] = sub_comlex(a[i], t);
            a[i] = add_comlex(a[i], t);
            w = mul_comlex(w, wn);
        }
    }
}

void fft1(Complex* a, int n, int inverse) {
    const int logn = log2(n);
    
    for (int i = 0; i < n; i++) {
        int k = 0;
        for (int l = 0; l < logn; l++)
            k |= ((i >> l) & 1) << (logn - l - 1);
        if (i < k) {
            Complex temp = a[i];
            a[i] = a[k];
            a[k] = temp;
        }
    }

    Complex wn = create_complex(cos(inverse * DPI / n), sin(inverse * DPI / n));    
    solve(a, n, wn);
}

PyObject* fft(PyObject* self, PyObject* args) {
    PyObject* list_obj;
    int inverse;

    if (!PyArg_ParseTuple(args, "Oi", &list_obj))
        return NULL;

    int n = PyList_Size(list_obj);
    Complex* a = (Complex*)malloc(n * sizeof(Complex));

    for (int i = 0; i < n; i++) {
        PyObject* item = PyList_GetItem(list_obj, i);
        double real = PyComplex_RealAsDouble(item);
        double imag = PyComplex_ImagAsDouble(item);
        a[i] = create_complex(real, imag);
    }

    fft1(a, n, -1);

    PyObject* result_list = PyList_New(n);
    for (int i = 0; i < n; i++) {
        PyObject* complex_obj = PyComplex_FromDoubles(a[i].real, a[i].imag);
        PyList_SetItem(result_list, i, complex_obj);
    }

    free(a);
    return result_list;
}

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

PyObject* hart(PyListObject* self, PyObject* args) {
    PyObject* list;

    if (!PyArg_ParseTuple(args, "O", &list))
        return NULL;

    if (!PyList_Check(list)) {
        PyErr_SetString(PyExc_TypeError, "Argument must be a list");
        return NULL;
    }

    Py_ssize_t N = PyList_GET_SIZE(list);
    PyObject* output = PyList_New(N);
    int* nist = (float_t*)malloc(N * sizeof(float_t));

    for (Py_ssize_t i = 0; i < N; i++) {
        nist[i] = PyFloat_AsDouble(PyList_GET_ITEM(list, i));
    }

    float_t temp;

    for (Py_ssize_t i = 0; i < N; i++) {
        temp = 0.0;

        double angle = DPI * i / N;
        for (Py_ssize_t j = 0; j < N; j++) {
            temp += nist[j] * (cos(angle * j) + sin(angle * j));
        }

        PyList_SET_ITEM(output, i, PyFloat_FromDouble(temp));
    }

    return output;
}

PyObject* ihart(PyObject* self, PyObject* args) {
    PyObject* list;

    if (!PyArg_ParseTuple(args, "O", &list))
        return NULL;

    if (!PyList_Check(list)) {
        PyErr_SetString(PyExc_TypeError, "Argument must be a list");
        return NULL;
    }

    Py_ssize_t N = PyList_GET_SIZE(list);
    PyObject* output = PyList_New(N);
    int* nist = (float_t*)malloc(N * sizeof(float_t));

    for (Py_ssize_t i = 0; i < N; i++) {
        nist[i] = PyFloat_AsDouble(PyList_GET_ITEM(list, i));
    }

    float_t temp;

    for (Py_ssize_t i = 0; i < N; i++) {
        temp = 0.0;

        double angle = DPI * i / N;
        for (Py_ssize_t j = 0; j < N; j++) {
            temp += nist[j] * (cos(angle * j) + sin(angle * j));
        }

        PyList_SET_ITEM(output, i, PyFloat_FromDouble(temp / N));
    }

    return output;
}

#endif




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

