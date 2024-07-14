#ifndef VECTOR_H
#define VECTOR_H

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include "val.h"

std::vector<double> vec_swap(const std::vector<double>& y);

class array {
private:
    std::vector<double> vec;
public:
    array() = default;

    array(std::vector<double> vec) : vec(std::move(vec)) {}

    array(int n) : vec(n) {}

    ~array() = default;

    int size() const {
        return static_cast<int>(this->vec.size()); // Исправлено предупреждение о возможной потере данных
    }

    void append(double val) {
        this->vec.push_back(val);
    }

    double front() const {
        return this->vec.front();
    }

    double back() {
        return this->vec.back();
    }

    double at(int index) {
        return this->vec.at(index);
    }

    array operator+ (array& vac) const {
        array vec_new;
        for (int i = 0; i < this->vec.size(); i++) {
            vec_new.append(this->vec[i]);
        }
        for (int i = 0; i < vac.vec.size(); i++) {
            vec_new.append(vac.vec[i]);
        }
        return vec_new;
    }

    array operator* (double& val) const {
        array vec_new;
        for (int i = 0; i < this->vec.size(); i++) {
            vec_new[i] = vec_new[i] * val;
        }
        return vec_new;
    }

    double& operator[](int index) {
        return this->vec[index];
    }

    const double& operator[](int index) const {
        return this->vec[index];
    }

    void reverse() {
        std::vector<double> temp;
        for (int i = this->vec.size() - 1; i >= 0; i--) {
            temp.push_back(this->vec[i]);
        }
        this->vec = temp;
    }

    friend std::ostream& operator<<(std::ostream& os, const array& arr);
};

#endif