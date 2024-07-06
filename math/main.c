#include "harmonic_c.h"

static PyObject* sum_list(PyObject* self, PyObject* args)
{
    Py_ssize_t i, n;
    long total = 0, value;
    PyObject* item;
    PyObject* list;

    if (!PyArg_ParseTuple(args, "O", &list))
        return NULL;

    n = PyList_Size(list);
    if (n < 0)
        return -1; /* Not a list */
    for (i = 0; i < n; i++) {
        item = PyList_GetItem(list, i); /* Can't fail */
        if (!PyLong_Check(item)) continue; /* Skip non-integers */
        value = PyLong_AsLong(item);
        if (value == -1 && PyErr_Occurred())
            /* Integer too big to fit in a C long, bail out */
            return -1;
        total += value;
    }
    return Py_BuildValue("i", total);
}

static PyMethodDef FtMethods[] = {
    {"ft", ft, METH_VARARGS, "Return the input list"},
    {"sum_list", sum_list, METH_VARARGS, "Sum the elements of a list"},
    {"fftn", fftn, METH_VARARGS, "Sum the elements of a list"},
    {"py_fft", py_fft, METH_VARARGS, "Sum the elements of a list"},
    {NULL, NULL, 0, NULL}
};

// Определение модуля
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "ftmodule",
    NULL,
    -1,
    FtMethods
};

PyMODINIT_FUNC PyInit_ftmodule(void) {
    return PyModule_Create(&moduledef);
}