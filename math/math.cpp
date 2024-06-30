#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "interpol.h"
#include "parser.h"
#include <pybind11/complex.h>
#include <cmath>
#include <stack>
#include "harmonic.h"
#include "vector.h"

namespace py = pybind11;

/*py::object find_index(const std::vector<py::object>& items, const py::object& value) {
    for (int i = 0; i < items.size(); i++) {
        if (items[i] == value) {
            return py::cast(i);
        }
    }   
    return py::none();
}*/


PYBIND11_MODULE(numera, m) {
    m.doc() = "";

    m.attr("pi") = py::float_(PI);

    // Define submodules
    py::module submodule = m.def_submodule("harmonic", "");
    submodule.def("fft", &fft, pybind11::arg("y"), "");
    submodule.def("ifft", &ifft, pybind11::arg("Y"), "");
    submodule.def("hart", &hart, pybind11::arg("y"), "");
    submodule.def("ihart", &ihart, pybind11::arg("H"), "");
    py::class_<array>(m, "array")
        .def(py::init<>())
        .def(py::init<std::vector<double>>())
        .def("size", &array::size)
        .def("append", &array::append)
        .def("front", &array::front)
        .def("back", &array::back)
        .def("at", &array::at)
        .def("reverse", &array::reverse)
        .def("__len__", &array::size)
        .def("__getitem__", [](const array& a, int i) {
        return a[i];
            })
        .def("__setitem__", [](array& a, int i, double v) {
                a[i] = v;
            })
                .def("__add__", &array::operator+)
                .def("__str__", [](const array& a) {
                std::ostringstream oss;
                oss << a;
                return oss.str();
                    })
                .def("__repr__", [](const array& a) {
                        std::ostringstream oss;
                        oss << a;
                        return oss.str();
                    });

    // Add functions to submodule
    //submodule.def("add", &add, "Addition function");
    //submodule.def("subtract", &subtract, "Subtraction function");
	//submodule.def("parser", &parser, "Subtraction function");
    //submodule.def("find_index", &find_index, "find_index function");
}