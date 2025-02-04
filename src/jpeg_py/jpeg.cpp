// Note it is recommended to include Python first, even before standard headers
#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <array>
#include <exception>
#include <stdexcept>
#include <string>

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

/// @brief Validates the quality factor is between 0 and 100 inclusive
/// @param qf Quality factor
void validate_quality_factor(int qf)
{
    if (qf<0)
    {
        throw_python_exception(PyExc_ValueError, "Quality factor cannot be less than 0");
    }
    else if (qf>100)
    {
        throw_python_exception(PyExc_ValueError, "Quality factor cannot be greater than 100");
    }
}

/// @brief Gets the subsampling used
/// @param ss_str Subsampling as a string
/// @return Subsampling as a JPEG::Subsampling
JPEG::Subsampling convert_subsampling_str(const std::string& ss_str)
{
    JPEG::Subsampling ss;

    if (ss_str=="4:4:4")
    {
        ss = JPEG::Subsampling::ss_4_4_4;
    }
    else if (ss_str=="4:2:0")
    {
        ss = JPEG::Subsampling::ss_4_2_0;
    }
    else
    {
        throw_python_exception(PyExc_ValueError, "Invalid subsampling");
    }

    return ss;
}

/// @brief Wraps the C++ JPEG:encode_greyscale
/// @param np_array Image data as a numpy array
/// @param qf Quality factor
/// @return Encoded image as a Python bytes object
PyObject* encode_greyscale_wrapper(PyObject* np_array, int qf)
{
    JPEG::Array_2d<double> array_2d;
    std::vector<unsigned char> encoded_image;

    array_2d = convert_numpy_to_array_2d(np_array);

    validate_quality_factor(qf);

    // Perform the encoding
    encoded_image = JPEG::encode_greyscale_image(array_2d, qf);

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

PyDoc_STRVAR(encode_greyscale_docstring,
    "encode_greyscale(data, qf=50)\n"   // Include function's signature first 
    "--\n\n"                            // We need this "--\n\n" so Python knows the first line is the function signature
    "\n"
    "Encodes a greyscale image as a JPEG\n"
    "\n"
    "Parameters\n"
    "----------\n"
    "data : numpy.ndarray\n"
    "    Image data to encode\n"
    "qf : int, optional\n"
    "    Quality factor from 0 to 100 inclusive used to set image quality using\n"
    "    the Independent JPEG Group's algorithm. A value of 50 corresponds to\n"
    "    using the JPEG specification's suggested quantization tables. The\n"
    "    default is 50.\n"
    "\n"
    "Returns\n"
    "-------\n"
    "bytes\n"
    "    The encoded image as a series of bytes"
);

static PyObject *encode_greyscale(PyObject *self, PyObject *args, PyObject* kwargs) {
    PyObject* np_array;
    int qf{ 50 };  // default value correpsonds to JPEG spec suggested tables

    // Names of all arguments, including positional and keyword
    const char* keywords[] = {
        "data",
        "qf",
        nullptr
    };

    // Note PyArg_ParseTupleAndKeywords() also ensures qf will fit in a C int. If not it raises a 
    // Python exception for up
    // Apparently using const_cast is appropriate here, PyArg_ParseTupleAndKeywords() shouldn't alter
    // keywords
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|i", const_cast<char**>(keywords), &np_array, &qf)) {
        return nullptr;
    }

    PyObject* encoded_image{};

    try
    {
        encoded_image = encode_greyscale_wrapper(np_array, qf);
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

/// @brief Wraps the C++ JPEG:encode_colour
/// @param red_np Red image component as a numpy array
/// @param green_np Green image component as a numpy array
/// @param blue_np Blue image component as a numpy array
/// @param qf Quality factor
/// @param ss_str Subsampling
/// @return Encoded image as a Python bytes object
PyObject* encode_colour_wrapper(
    PyObject* red_np, PyObject* green_np, PyObject* blue_np, int qf, const std::string& ss_str
)
{
    JPEG::Array_2d<double> red, green, blue;

    red = convert_numpy_to_array_2d(red_np);
    green = convert_numpy_to_array_2d(green_np);
    blue = convert_numpy_to_array_2d(blue_np);

    // Verify components have the same shape
    if (
         red.shape()[0]!=green.shape()[0] ||  red.shape()[1]!=green.shape()[1] ||
        blue.shape()[0]!=green.shape()[0] || blue.shape()[1]!=green.shape()[1]
        )
    {
        throw_python_exception(PyExc_ValueError, "Components did not all have the same shape");
    }

    validate_quality_factor(qf);

    JPEG::Subsampling ss{ convert_subsampling_str(ss_str) };

    // Perform the encoding
    std::vector<unsigned char> encoded_image;

    encoded_image = JPEG::encode_colour_image(red, green, blue, qf, ss);

    // Now convert to a bytes object
    // note PyBytes_FromStringAndSize() returns null on failure and sets the Python error state
    PyObject* bytes{};

    bytes = PyBytes_FromStringAndSize((char*)encoded_image.data(), encoded_image.size());

    if (!bytes)
    {
        throw Python_Exception("Error creating bytes object");
    }

    return bytes;
}

PyDoc_STRVAR(encode_colour_docstring,
    // Include function's signature first 
    "encode_colour(red, green, blue, qf=50, ss='4:4:4')\n"   
    "--\n\n"                                     // We need this "--\n\n" so Python knows the first line is the function signature
    "\n"
    "Encodes a colour image as a JPEG\n"
    "\n"
    "Parameters\n"
    "----------\n"
    "red : numpy.ndarray\n"
    "    Red image component\n"
    "green : numpy.ndarray\n"
    "    Green image component\n"
    "blue : numpy.ndarray\n"
    "    Blue image component\n"
    "qf : int, optional\n"
    "    Quality factor from 0 to 100 inclusive used to set image quality using\n"
    "    the Independent JPEG Group's algorithm. A value of 50 corresponds to\n"
    "    using the JPEG specification's suggested quantization tables. The\n"
    "    default is 50.\n"
    "ss : str, optional\n"
    "    Subsampling to perform, either '4:4:4' (no subsampling) or '4:2:0'. The\n"
    "    default is '4:4:4'.\n"
    "\n"
    "Returns\n"
    "-------\n"
    "bytes\n"
    "    The encoded image as a series of bytes"
);

static PyObject *encode_colour(PyObject *self, PyObject *args, PyObject* kwargs) {
    PyObject* red;
    PyObject* green;
    PyObject* blue;
    int qf{ 50 };   // default value correpsonds to JPEG spec suggested tables
    const char* ss{ "4:4:4" };     // subsampling

    // Names of all arguments, including positional and keyword
    const char* keywords[] = {
        "red",
        "green",
        "blue",
        "qf",
        "ss",
        nullptr
    };

    // Note PyArg_ParseTupleAndKeywords() also ensures qf will fit in a C int. If not it raises a 
    // Python exception for up
    // Apparently using const_cast is appropriate here, PyArg_ParseTupleAndKeywords() shouldn't alter
    // keywords
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOO|is", const_cast<char**>(keywords), &red, &green, &blue, &qf, &ss)) {
        return nullptr;
    }

    // Set the subsampling string
    std::string ss_str{ ss };

    PyObject* encoded_image{};

    try
    {
        encoded_image = encode_colour_wrapper(red, green, blue, qf, ss_str);
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
    {"encode_greyscale", (PyCFunction)encode_greyscale, METH_VARARGS | METH_KEYWORDS, encode_greyscale_docstring},
    {"encode_colour", (PyCFunction)encode_colour, METH_VARARGS | METH_KEYWORDS, encode_colour_docstring},
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
