// Note it is recommended to include Python first, even before standard headers
#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <array>
#include <exception>
#include <stdexcept>

// Numpy array API
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/ndarrayobject.h"

#include "jpeg_cpp/array.h"
#include "jpeg_cpp/jpeg.h"

/// @brief Exception class to indicate a Python exception has been raised
/// and the Python error state needs to be checked
class Python_Exception : public std::exception
{
    const char* message;

public:
    /// @brief Constructs a Python_Exception instance
    /// @param message Error message, does not take ownership
    Python_Exception(const char* message) noexcept : message{message} {}

    const char* what() const noexcept override { return message; }
};

/// @brief Throes a Python_exception and sets the Python error state
/// @param exception Python exception to set
/// @param message Error message
void throw_python_exception(PyObject* exception, const char* message)
{
    PyErr_SetString(exception, message);
    throw Python_Exception(message);
}

/// @brief Constructs an Arrat_2d<double> instance from a NumPy array. Works by iterating 
/// using NumPy iteraotrs so should work regardless of any array slicing, array order 
/// (C/Fortran) etc. convert_numpy_to_array_2d() should be prefered.
/// @param np_array A 2d NumPy array with type NPY_UINT8
/// @return The newly constructed Arrat_2d<double> 
JPEG::Array_2d<double> convert_numpy_generic_to_array_2d(PyArrayObject* np_array)
{
    std::array<size_t, 2> shape{
        static_cast<size_t>(PyArray_DIMS(np_array)[0]),
        static_cast<size_t>(PyArray_DIMS(np_array)[1]),
    };

    JPEG::Array_2d<double> array_2d{ shape[0], shape[1] };

    NpyIter *iter{};
    NpyIter_IterNextFunc *iternext{};
    npy_intp multi_index[2];
    char** dataptr{};

    // Unfourtunately NpyIter_Deallocate() can itself raise an exception so we
    // can't wrap iter up as a smart pointer. Best we can do
    iter = NpyIter_New(
        np_array, NPY_ITER_READONLY | NPY_ITER_MULTI_INDEX | NPY_ITER_REFS_OK,
        NPY_KEEPORDER, NPY_NO_CASTING, NULL
    );

    if (iter==nullptr) {
        // NpyIter_New() sets Python error state, just need to return
        throw Python_Exception("Error creating NumPy iterator");
    }

    if (NpyIter_GetIterSize(iter) != 0) {
        iternext = NpyIter_GetIterNext(iter, NULL);
        dataptr = NpyIter_GetDataPtrArray(iter);

        if (iternext == NULL) {
            NpyIter_Deallocate(iter);
            throw Python_Exception("Error getting iternext");
        }

        NpyIter_GetMultiIndexFunc *get_multi_index = NpyIter_GetGetMultiIndex(iter, NULL);

        if (get_multi_index == NULL) {
            NpyIter_Deallocate(iter);
            throw Python_Exception("Error getting get_multi_index");
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
        throw Python_Exception("Error in NpyIter_Deallocate()");
    }

    return array_2d;
}

/// @brief Constructs an Arrat_2d<double> instance from a NumPy array.
/// @param obj NumPy array as a Python object
/// @return The newly constructed Arrat_2d<double> 
JPEG::Array_2d<double> convert_numpy_to_array_2d(PyObject* obj)
{
    // Check it is a numpy array
    // This allows for subclasses of numpy arrays, probably OK?
    if (!PyArray_Check(obj))
    {   
        throw_python_exception(PyExc_ValueError, "Argument was not a NumPy array");
    }

    PyArrayObject* np_array{ (PyArrayObject*)obj };

    // Check it has correct type
    if (PyArray_TYPE(np_array)!=NPY_UINT8)
    {
        throw_python_exception(PyExc_ValueError, "Array did not have expected type, expected np.uint8");
    }

    // Check it is 2d
    if (PyArray_NDIM(np_array)!=2)
    {
        throw_python_exception(PyExc_ValueError, "Array did not have expected number of dimensions, expected a 2d array");
    }

    // Check its size is greater than zero in both directions
    std::array<size_t, 2> shape{
        static_cast<size_t>(PyArray_DIMS(np_array)[0]),
        static_cast<size_t>(PyArray_DIMS(np_array)[1]),
    };

    if (shape[0]==0 || shape[1]==0)
    {
        throw_python_exception(PyExc_ValueError, "Array cannot have size zero in any dimension");
    }

    // This works for a generic numpy array, regardless of how it exists in memory
    // In future we might be able to get a performance improvement for arrays that
    // are C contiguous. We can use PyArray_IS_C_CONTIGUOUS() to check this. Note 
    // this question regarding strides: https://github.com/numpy/numpy/issues/15979
    return convert_numpy_generic_to_array_2d(np_array);
}

PyDoc_STRVAR(encode_greyscale_docstring,
    "encode_greyscale(data)\n"  // Include function's signature first 
    "--\n\n"                    // We need this "--\n\n" so Python knows the first line is the function signature
    "\n"
    "Encodes a greyscale image as a JPEG\n"
    "\n"
    "Parameters\n"
    "----------\n"
    "data : numpy.ndarray\n"
    "    Image data to encode\n"
    "\n"
    "Returns\n"
    "-------\n"
    "bytes\n"
    "    The encoded image as a series of bytes"
);

/// @brief Wraps the C++ JPEG:encode_greyscale
/// @param np_array Image data as a numpy array
/// @return Encoded image as a Python bytes object
PyObject* encode_greyscale_wrapper(PyObject* np_array)
{
    JPEG::Array_2d<double> array_2d;
    std::vector<unsigned char> encoded_image;

    array_2d = convert_numpy_to_array_2d(np_array);

    // Perform the encoding
    encoded_image = JPEG::encode_greyscale_image(array_2d);

    // Now convert to a bytes object
    // note PyBytes_FromStringAndSize() returns null on failure
    PyObject* bytes{};

    bytes = PyBytes_FromStringAndSize((char*)encoded_image.data(), encoded_image.size());

    if (!bytes)
    {
        throw Python_Exception("Error creating bytes object");
    }

    return bytes;
}

static PyObject *encode_greyscale(PyObject *self, PyObject *args) {
    PyObject* np_array;

    if (!PyArg_ParseTuple(args, "O", &np_array)) {
        return nullptr;
    }

    PyObject* encoded_image{};

    try
    {
        encoded_image = encode_greyscale_wrapper(np_array);
    }
    catch(const Python_Exception& e)
    {
        // Something raised a Python exception. Python exception state is 
        // already set so just return
        return nullptr;
    }
    catch(const std::exception& e)
    {
        // Unexpected error, Python exception state has not been set
        // Set the Python exception state as a runtime error
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return nullptr;
    }
    
    return encoded_image;
}

static PyMethodDef jpeg_methods[] = {
    {"encode_greyscale", encode_greyscale, METH_VARARGS, encode_greyscale_docstring},
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
