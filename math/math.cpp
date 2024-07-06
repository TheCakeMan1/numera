#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "interpol.h"
#include "parser.h"
#include <pybind11/complex.h>
#include <cmath>
#include <stack>
#include "harmonic.h"
#include "vector.h"
#include "ni.h"
#include "nf.h"
#include "nui.h"
extern "C" {
#include "vector_c.h"
}

namespace py = pybind11;


/*py::object find_index(const std::vector<py::object>& items, const py::object& value) {
    for (int i = 0; i < items.size(); i++) {
        if (items[i] == value) {
            return py::cast(i);
        }
    }   
    return py::none();
}*/

int add(int a, int b) {
    return a + b;
}

PYBIND11_MODULE(numera, m) {
    m.doc() = "";

    m.attr("pi") = py::float_(PI);

    // Define submodules
    py::module submodule = m.def_submodule("harmonic", "");
    submodule.def("add", &add, "");
    submodule.def("fft", &fft, pybind11::arg("y"), "");
    submodule.def("ift", &ift, pybind11::arg("Y"), "");
    submodule.def("hart", &hart, pybind11::arg("y"), "");
    submodule.def("ihart", &ihart, pybind11::arg("H"), "");
    submodule.def("convolution", &convolution, pybind11::arg("y_1"), pybind11::arg("y_2"), "");


    py::class_<array>(m, "array")
        .def(py::init<>())
        .def(py::init<std::vector<double>>())
        .def("size", &array::size, R"pbdoc(
This function will allow you to find out the dimension of the vector.
            
Parameters
==========
void
            
Returns
=======
Returns an integer indicating the length of the array.

Examples
========
>>> import numera as nu
>>> vector = array([1, 2, 3])
>>> print(vector)
        )pbdoc")
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
    py::class_<ni>(m, "ni")
        .def(py::init<>())
        .def(py::init<int>())
        .def(py::init<float>())
        .def(py::init<double>())
        .def(py::init<const ni&>())
        .def("toInt", &ni::toInt)
        .def("toFloat", &ni::toFloat)
        .def("__add__", static_cast<ni(ni::*)(const ni&) const>(&ni::operator+), py::is_operator())
        .def("__add__", static_cast<ni(ni::*)(int) const>(&ni::operator+), py::is_operator())
        .def("__add__", static_cast<ni(ni::*)(double) const>(&ni::operator+), py::is_operator())
        //.def("__add__", static_cast<ni(ni::*)(const nf&) const>(&ni::operator+), py::is_operator())
        .def("__sub__", &ni::operator-, py::is_operator())
        .def("__mul__", &ni::operator*, py::is_operator())
        .def("__str__", &ni::print)
        .def("__repr__", &ni::print);

    py::class_<nui>(m, "nui")
        .def(py::init<>())
        .def(py::init<int>())
        .def(py::init<float>())
        .def(py::init<double>())
        .def(py::init<const nui&>())
        .def("toInt", &nui::toInt)
        .def("toFloat", &nui::toFloat)
        .def("__add__", static_cast<nui(nui::*)(const nui&) const>(&nui::operator+), py::is_operator())
        .def("__add__", static_cast<nui(nui::*)(int) const>(&nui::operator+), py::is_operator())
        .def("__add__", static_cast<nui(nui::*)(double) const>(&nui::operator+), py::is_operator())
        //.def("__add__", static_cast<ni(ni::*)(const nf&) const>(&ni::operator+), py::is_operator())
        .def("__sub__", &nui::operator-, py::is_operator())
        .def("__mul__", &nui::operator*, py::is_operator())
        .def("__str__", &nui::print)
        .def("__repr__", &nui::print);

    py::class_<nf>(m, "nf")
        .def(py::init<>())
        .def(py::init<float>())
        //.def(py::init<ni>())
        .def(py::init<const nf&>())
        .def("toFloat", &nf::toFloat);

    py::class_<bint>(m, "Bint")
        .def(py::init<>())
        .def_readwrite("array", &bint::array)
        .def_readwrite("size", &bint::size)
        .def_readwrite("size_last", &bint::size_last)
        .def_readwrite("capacity", &bint::capacity)
        .def_readwrite("negative", &bint::negative);

    m.def("initVector", &initVector, "Initialize bint from string");
    m.def("print_bint", &print_bint, "Print bint object");
    m.def("add_bint", &add_bint, "Add two bint objects");
        

    // Add functions to submodule
    //submodule.def("add", &add, "Addition function");
    //submodule.def("subtract", &subtract, "Subtraction function");
	//submodule.def("parser", &parser, "Subtraction function");
    //submodule.def("find_index", &find_index, "find_index function");
}