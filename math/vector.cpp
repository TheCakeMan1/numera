#include "vector.h"

std::ostream& operator<<(std::ostream& os, const array& arr) {
    os << "[";
    for (size_t i = 0; i < arr.vec.size(); ++i) {
        os << arr.vec[i];
        if (i != arr.vec.size() - 1) {
            os << ", ";
        }
    }
    os << "]";
    return os;
}

std::vector<double> vec_swap(const std::vector<double>& y) {
    std::vector<double> temp;
    for (auto i : y) {
        temp.push_back(i);
    }
    return temp;
}

array zeros(int n){
    return array(std::vector<double>(n, 0.0));
}