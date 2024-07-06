#include <Python.h>

// Эту функцию мы вызываем из Python кода
static PyObject* addList_add(PyObject* self, PyObject* args) {
    PyObject* listObj;
    if (!PyArg_ParseTuple(args, "O", &listObj))
        return NULL;

    long length = PyList_Size(listObj);

    // Проходимся по всем элементам
    Py_ssize_t i;
    PyObject* sum = PyLong_FromLong(0);
    for (i = 0; i < length; i++) {
        PyObject* temp = PyList_GetItem(listObj, i);

        PyObject* elem = PyLong_AsLong(temp);
        sum = PyNumber_Add(sum, temp);
    }

    // Возвращаемое в Python-код значение также Python-объект
    // Приводим C long к Python integer
    return Py_BuildValue("i", sum);
}

// Немного документации для `add`
static char addList_docs[] =
"add( ): add all elements of the list\n";

/*
Эта таблица содержит необходимую информацию о функциях модуля
<имя функции в модуле Python>, <фактическая функция>,
<ожидаемые типы аргументов функции>, <документация функции>
*/
static PyMethodDef addList_funcs[] = {
    {"add", (PyCFunction)addList_add, METH_VARARGS, addList_docs},
    {NULL, NULL, 0, NULL}
};

/*
addList имя модуля и это блок его инициализации.
<желаемое имя модуля>, <таблица информации>, <документация модуля>
*/

static struct PyModuleDef addList_module = {
    PyModuleDef_HEAD_INIT,
    "addList",
    "Add all ze lists",
    -1,
    addList_funcs
};

PyMODINIT_FUNC PyInit_addList(void) {
    return PyModule_Create(&addList_module);
}
