#ifndef ENCODING_H
#define ENCODING_H

#include <array>
#include <cassert>

// Needed to access pi constant
#define _USE_MATH_DEFINES
#include <cmath>

#include <stdexcept>
#include <vector>

#include "jpeg_cpp/array.h"
#include "jpeg_cpp/bit_string.h"
#include "jpeg_cpp/huff_table.h"
#include "jpeg_cpp/q_table.h"

namespace JPEG
{
    /// @brief Indices for zig-zag ordering for a data unit stored in row-major order
    constexpr std::array<size_t, 64> zig_zag_order{ 
        0, 1, 8, 16, 9, 2, 3, 10, 17, 24, 32, 
        25, 18, 11, 4, 5, 12, 19, 26, 33, 40, 
        48, 41, 34, 27, 20, 13, 6, 7, 14, 21, 
        28, 35, 42, 49, 56, 57, 50, 43, 36, 
        29, 22, 15, 23, 30, 37, 44, 51, 58, 
        59, 52, 45, 38, 31, 39, 46, 53, 60, 
        61, 54, 47, 55, 62, 63
    };

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

    /// @brief Computes ssss value of a number. For n!=0, this is the number such that 
    /// n >> ssss == 1. For n==0 this is 0. This function assumes ssss<=11, i.e. n<=2047.
    /// @param n 
    /// @return The ssss of n
    unsigned int compute_ssss(unsigned int n);

    /// @brief Encodes a DC coefficient and appends the result to the Bit_String
    /// @param bs Bit_String to append encoded DC coefficient to
    /// @param diff Difference between the current DC coefficient and the previous DC coefficient
    /// @param huff_table DC Huffman table
    void encode_DC_coeff(Bit_String& bs, int diff, const Huff_Table& huff_table);

    /// @brief Encodes a non-zero coefficient
    /// @param bs Bit_String to append encoded AC coefficient to
    /// @param coeff Current AC coefficient to encode
    /// @param rrrr Current run length of zeros
    /// @param huff_table AC Huffman table
    void encode_AC_coeff(Bit_String& bs, int coeff, unsigned int rrrr, const Huff_Table& huff_table);

    /// @brief Encodes the AC coefficients of a data unit and appends the result to the Bit_String
    /// @param bs Bit_String to append encoded DC coefficient to
    /// @param du_array DU_Array with the data unit to encode
    /// @param du_ind Index of the data unit within the DU_Array
    /// @param huff_table AC Huffman table
    void encode_AC_coeffs(Bit_String& bs, const DU_Array<double>& du_array, size_t du_ind,
                            const Huff_Table& huff_table);

    /// @brief Encodes a data unit using the sequential mode
    /// @param bs Bit_String to append encoded DC coefficient to
    /// @param du_array DU_Array with the data unit to encode
    /// @param du_ind Index of the data unit within the DU_Array
    /// @param prev_dc DC coefficient of the previously encoded data unit
    /// @param huff_table_dc DC Huffman table
    /// @param huff_table_ac AC Huffman table
    void encode_data_unit_sequential(Bit_String& bs, const DU_Array<double>& du_array, size_t du_ind, int prev_dc,
                        const Huff_Table& huff_table_dc, const Huff_Table& huff_table_ac);

    /// @brief Appends a marker segment describing the given quantization tables
    /// @param out Output to append encoded marker segment to
    /// @param q_tables Quantization tables to encode
    /// @param destination_indices Destination indices of the tables in q_tables
    void append_q_table_marker_segment(std::vector<unsigned char>& out, const std::vector<Q_Table>& q_tables, 
                        std::vector<unsigned int>& destination_indices);
}

#endif
