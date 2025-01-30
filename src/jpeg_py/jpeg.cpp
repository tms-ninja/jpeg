// Note it is recommended to include Python first, even before standard headers
#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <array>
#include <stdexcept>

// Numpy array API
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/ndarrayobject.h"

#include "jpeg_cpp/array.h"
#include "jpeg_cpp/jpeg.h"

float square(float x) { return x * x; }

static PyObject *square_wrapper(PyObject *self, PyObject *args) {
    float input, result;
    if (!PyArg_ParseTuple(args, "f", &input)) {
    return NULL;
    }
    result = square(input);
    return PyFloat_FromDouble(result);
}

int convert_numpy_to_array_2d(PyObject* obj, JPEG::Array_2d<double>& array_2d)
{
    // Check it is a numpy array
    // This allows for subclasses of numpy arrays, probably OK?
    if (!PyArray_Check(obj))
    {
        PyErr_SetString(PyExc_ValueError, "Argument was not a NumPy array");
        return -1;
    }

    PyArrayObject* np_array{ (PyArrayObject*)obj };

    // Check it has correct type
    if (PyArray_TYPE(np_array)!=NPY_UINT8)
    {
        PyErr_SetString(PyExc_ValueError, "Array did not have expected type, expected np.uint8");
        return -1;
    }

    // Check it is 2d
    if (PyArray_NDIM(np_array)!=2)
    {
        PyErr_SetString(PyExc_ValueError, "Array did not have expected number of dimensions, expected a 2d array");
        return -1;
    }

    // Check its size is greater than zero in both directions
    std::array<size_t, 2> shape{
        static_cast<size_t>(PyArray_DIMS(np_array)[0]),
        static_cast<size_t>(PyArray_DIMS(np_array)[1]),
    };

    if (shape[0]==0 || shape[1]==0)
    {
        PyErr_SetString(PyExc_ValueError, "Array cannot have size zero in any dimension");
        return -1;
    }

    // Should be able to copy it
    array_2d = JPEG::Array_2d<double>{ shape[0], shape[1] };

    // This all uses 
    NpyIter *iter{};
    NpyIter_IterNextFunc *iternext{};
    npy_intp multi_index[2];
    char** dataptr{};

    iter = NpyIter_New(
        np_array, NPY_ITER_READONLY | NPY_ITER_MULTI_INDEX | NPY_ITER_REFS_OK,
        NPY_KEEPORDER, NPY_NO_CASTING, NULL
    );

    if (iter==nullptr) {
        // NpyIter_New() sets Python error state, just need to return
        return -1;
    }

    if (NpyIter_GetIterSize(iter) != 0) {
        iternext = NpyIter_GetIterNext(iter, NULL);
        dataptr = NpyIter_GetDataPtrArray(iter);

        if (iternext == NULL) {
            NpyIter_Deallocate(iter);
            return -1;
        }

        NpyIter_GetMultiIndexFunc *get_multi_index = NpyIter_GetGetMultiIndex(iter, NULL);

        if (get_multi_index == NULL) {
            NpyIter_Deallocate(iter);
            return -1;
        }

        do {
            // Gets the (i, j) index of the next element in the numpy array
            get_multi_index(iter, multi_index);
            
            // Pointer to the element in the nump array
            char* data{ *dataptr };

            array_2d(multi_index[0], multi_index[1]) = static_cast<double>(*(unsigned char*)data);

        } while (iternext(iter));
    }

    if (!NpyIter_Deallocate(iter)) {
        return -1;
    }

    return 0;
}

static PyObject *encode_greyscale(PyObject *self, PyObject *args) {
    PyObject* np_array;

    if (!PyArg_ParseTuple(args, "O", &np_array)) {
        return nullptr;
    }

    // Initalize to something temporary for now
    JPEG::Array_2d<double> array_2d{ 1, 1 };

    if (convert_numpy_to_array_2d(np_array, array_2d))
    {
        // Return null pointer to indicate an error to python
        return nullptr;
    }

    // Perform the encoding
    std::vector<unsigned char> encoded_image{ JPEG::encode_greyscale_image(array_2d) };

    // Now convert to a bytes object
    // note PyBytes_FromStringAndSize() returns null on failiure so can just return
    return PyBytes_FromStringAndSize((char*)encoded_image.data(), encoded_image.size());
}

static PyMethodDef jpeg_methods[] = {
    {"square", square_wrapper, METH_VARARGS, "Square function"},
    {"encode_greyscale", encode_greyscale, METH_VARARGS, "Encodes a greyscale image"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef jpeg_module = {
    PyModuleDef_HEAD_INIT, "jpeg",
    NULL, -1, jpeg_methods
};

/* name here must match extension name, with PyInit_ prefix */
PyMODINIT_FUNC PyInit_jpeg(void) {
    // Initalize numpy C API
    if (PyArray_ImportNumPyAPI() < 0) {
    return nullptr;
    }

    return PyModule_Create(&jpeg_module);
}
