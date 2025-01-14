#ifndef ENCODING_H
#define ENCODING_H

#include <array>

// Needed to access pi constant
#define _USE_MATH_DEFINES
#include <cmath>

#include "jpeg_cpp/array.h"
#include "jpeg_cpp/q_table.h"

namespace JPEG
{
    void apply_level_shift(DU_Array<double>& array);

    /// @brief Loads matrix for computing the discrete cosine transform
    /// This is defined as M_ij = 1/2 cos((2*i+1)*j*pi/16)
    /// @return DCT matrix
    Array_2d<double> load_DCT_matrix();

    /// @brief Loads the transpose of the DCT matrix
    /// @return Transpose of the DCT matrix
    Array_2d<double> load_DCT_matrix_transpose();

    const static Array_2d<double> DCT_matrix{ load_DCT_matrix() };
    const static Array_2d<double> DCT_matrix_transpose{ load_DCT_matrix_transpose() };

    /// @brief Computes the matrix multiplication c = a @ b. It assumes matrices are square
    /// and are stored in row-major order.
    /// @param a 
    /// @param b 
    /// @param c 
    /// @param N Size of the matrix, ie. a, b c have shape (N, N)
    void mat_mul(const double* a, const double* b, double* c, size_t N);

    void apply_DCT(DU_Array<double>& array);

    void apply_quantization(DU_Array<double>& array, const Q_Table& q_table);
}

#endif
