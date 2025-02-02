#ifndef ENCODING_H
#define ENCODING_H

// Needed to access pi constant
// Move to the top so MSVC can find M_PI
#define _USE_MATH_DEFINES
#include <cmath>

#include <array>
#include <cassert>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "jpeg_cpp/array.h"
#include "jpeg_cpp/bit_string.h"
#include "jpeg_cpp/component.h"
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

    /// @brief Appends the last two bytes to the vector, most significant first
    /// @param out Vector to append bytes to
    /// @param bytes Number containing bytes
    void append_two_bytes(std::vector<unsigned char>& out, unsigned long long bytes);

    /// @brief Appends a composite byte to the vector
    /// @param out Vector to append byte to
    /// @param high Source of the 4 high bits, taken from the 4 low bits of high, most significant first
    /// @param low Source of the 4 low bits, taken from the 4 low bits of low, most significant first
    void append_composite_byte(std::vector<unsigned char>& out, unsigned long long high, unsigned long long low);

    enum class Marker {
        Start_Of_Image = 0xFFD8,
        End_Of_Image = 0xFFD9,
        Define_Quantization_Table = 0xFFDB,
        Define_Huffman_Table = 0xFFC4,
        Start_Of_Frame_0_Baseline_DCT = 0xFFC0,
        Start_Of_Scan = 0xFFDA,
    };

    /// @brief Appends the marker to the output
    /// @param out Vector to append byte to
    /// @param marker JPEG marker
    void append_marker(std::vector<unsigned char>& out, Marker marker);

    /// @brief Applies the level shift of subtracting 128 to the array
    /// @param array Array to apply level shift to
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

    /// @brief Applies a quasi discrete cosine transform to the array. Specifically, it applies the transform
    /// less the factors of 1/sqrt(2) that apply to elements in the top row/left column
    /// @param array Array to apply the DCT to
    void apply_DCT(DU_Array<double>& array);

    /// @brief Quantizes the elements in the given array. It internally accounts for factors of 1/sqrt(2) that
    /// apply to elements in the top row/left column from the DCT stage. These do not need to be accounted for 
    /// by the caller.
    /// @param array Array to apply the DCT to
    /// @param q_table Quantization table used to quantize the elements in array
    void apply_quantization(DU_Array<double>& array, const Q_Table& q_table);

    /// @brief Computes ssss value of a number. For n!=0, this is the number such that 
    /// n >> ssss == 1. For n==0 this is 0. This function assumes ssss<=11, i.e. n<=2047.
    /// @param n Number to compute ssss of
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

    /// @brief Encodes a data unit using the sequential mode. Note this only performs entropy encoding.
    /// @param bs Bit_String to append encoded data unit to
    /// @param du_array DU_Array with the data unit to encode
    /// @param du_ind Index of the data unit within the DU_Array
    /// @param prev_dc DC coefficient of the previously encoded data unit
    /// @param huff_table_dc DC Huffman table
    /// @param huff_table_ac AC Huffman table
    /// @return DC coefficient of the current data unit
    int encode_data_unit_sequential(Bit_String& bs, const DU_Array<double>& du_array, size_t du_ind, int prev_dc,
                        const Huff_Table& huff_table_dc, const Huff_Table& huff_table_ac);

    /// @brief Appends the encoded MCU to the Bit_String. Note this only performs entropy encoding.
    /// @param bs Bit_String to append encoded MCU to
    /// @param prev_dc Previous DC coefficients of each component
    /// @param du_ind Index of the next data unit to encode for each component
    /// @param arrays Arrays comtaining the data units for each component
    /// @param comps Component descriptors for each component
    /// @param dc_tables List of DC tables. Note which table is used for each component is taken from comps
    /// @param ac_tables List of AC tables. Note which table is used for each component is taken from comps
    void append_mcu(Bit_String& bs, std::vector<int>& prev_dc, std::vector<size_t>& du_ind, const std::vector<DU_Array<double>>& arrays, 
        const std::vector<Comp_Info>& comp_infos, const std::vector<Huff_Table>& dc_tables, const std::vector<Huff_Table>& ac_tables);

    /// @brief Appends a marker segment describing the given quantization tables
    /// @param out Output to append encoded marker segment to
    /// @param q_tables Quantization tables to encode
    /// @param destination_indices Destination indices of the tables in q_tables
    void append_q_table_marker_segment(std::vector<unsigned char>& out, const std::vector<Q_Table>& q_tables, 
                        std::vector<unsigned int>& destination_indices);
    
    /// @brief Represents a reference to a Huffman table along with information like its DC/AC type
    struct Huff_Table_Ref
    {
        enum class Huff_Table_Type
        {
            AC,
            DC
        };

        const Huff_Table& table;
        Huff_Table_Type type;
        size_t destination_ind;

        Huff_Table_Ref(const Huff_Table& huff_table, Huff_Table_Type table_type, size_t destination) 
        : table{ huff_table }, type{ table_type }, destination_ind{ destination }
        {

        }
    };

    /// @brief Appends the BITS list giving the number of codes of lengths 1 to 16
    /// @param out Output to append encoded BITS list to
    /// @param huff_table Table to encode
    void append_BITS_array(std::vector<unsigned char>& out, const Huff_Table& huff_table);

    /// @brief Appends the HUFFVAL list giving the symbols in the Huffman table. Note this 
    /// assumes the Huffman table uses cannonical Huffman encoding. That, is for any i<j and 
    /// table[i].size()==table[i].size() then table[i]<table[j].
    /// @param out Output to append encoded HUFFVAL list to
    /// @param huff_table Table to encode
    void append_HUFFVAL_array(std::vector<unsigned char>& out, const Huff_Table& huff_table);

    /// @brief Appends encoding for a specific Huffman table to the output. This only includes
    /// the data for this specific Huffman table, that is the table type, destination and the 
    /// BITS/HUFFVAL arrays. It does not include the DHT marker or Huffman table definition 
    /// length.
    /// 
    /// Note this assumes the Huffman table uses cannonical Huffman encoding. That is, for any 
    /// i<j and table[i].size()==table[i].size() then table[i]<table[j].
    /// @param out Output to append encoded Huffman table data to
    /// @param huff_table Table to encode
    void append_huff_table_data(std::vector<unsigned char>& out, const Huff_Table_Ref& huff_table);

    /// @brief Appends Huffman table specification marker segment to the output. Note this 
    /// assumes the Huffman table uses cannonical Huffman encoding. That is, for any i<j 
    /// and table[i].size()==table[i].size() then table[i]<table[j].
    /// @param out Output to append encoded Huffman table data to
    /// @param tables Tables to encode
    void append_huff_table_marker_segment(std::vector<unsigned char>& out, std::vector<Huff_Table_Ref> tables);

    /// @brief Appends a scan header to the output.
    /// @param out Output to append scan header to
    /// @param comp_infos List of components in the scan
    void append_scan_header(std::vector<unsigned char>& out, const std::vector<Comp_Info>& comp_infos);

    /// @brief Appends a frame header to the output.
    /// @param out Output to append frame header to
    /// @param Y Height of image
    /// @param X Width of image
    /// @param comp_infos List of components in the frame
    void append_frame_header(std::vector<unsigned char>& out, unsigned int Y, unsigned int X, 
        const std::vector<Comp_Info>& comp_infos);

    /// @brief Encodes a scan. Note does not perform level shift/DCT/quantization
    /// @param out Output to append scan to
    /// @param arrays Arrays comtaining the data units for each component. Level shift/DCT/quantization should have already
    /// been performed
    /// @param comp_infos List of components in the scan
    /// @param dc_tables List of DC tables. Note which table is used for each component is taken from comp_infos
    /// @param ac_tables List of AC tables. Note which table is used for each component is taken from comp_infos
    /// @param q_tables List of quantization tables. Note which table is used for each component is taken from comp_infos
    void encode_scan(std::vector<unsigned char>& out, std::vector<DU_Array<double>>& arrays, 
        const std::vector<Comp_Info>& comp_infos, const std::vector<Huff_Table>& dc_tables, const std::vector<Huff_Table>& ac_tables);

    /// @brief Encodes a frame. Note does not perform level shift/DCT/quantization
    /// @param out Output to append frame to
    /// @param Y Height of the image in pixels
    /// @param X Width of the image in pixels
    /// @param arrays Arrays comtaining the data units for each component. Level shift/DCT/quantization should have already
    /// been performed
    /// @param comp_infos List of components in the scan
    /// @param dc_tables List of DC tables. Note which table is used for each component is taken from comp_infos
    /// @param ac_tables List of AC tables. Note which table is used for each component is taken from comp_infos
    /// @param q_tables List of quantization tables. Note which table is used for each component is taken from comp_infos
    void encode_frame(std::vector<unsigned char>& out, unsigned int Y, unsigned int X, std::vector<DU_Array<double>>& arrays, 
        const std::vector<Comp_Info>& comp_infos, const std::vector<Huff_Table>& dc_tables, const std::vector<Huff_Table>& ac_tables,
        const std::vector<Q_Table>& q_tables);

    /// @brief Converts from a 2d array representation to a data unit array representation. Note the array_2d should be 
    /// appropriately enlarged to ensure an integer number of MCUs
    /// @param array_2d Array_2d to convert
    /// @param comp_info Comp_Info describing the component, including the H and V sampling factors
    /// @return The Array_2d in DU_Array representation
    DU_Array<double> convert_to_DU_Array(const Array_2d<double>& array_2d, const Comp_Info& comp_info);

    /// @brief Encodes an image, note components should be appropriately enlarged to ensure an integer number of MCUs
    /// @param Y Height of the image in pixels
    /// @param X Width of the image in pixels
    /// @param arrays Image components. Any colour transformation/subsampling should have already be performed. In particular,
    /// components should be appropriately enlarged to ensure an integer number of MCUs
    /// @param comp_infos Comp_Infos for each corresponding component
    /// @param dc_tables DC Huffman tables to be used in the encoding process
    /// @param ac_tables AC Huffman tables to be used in the encoding process
    /// @param q_tables Quantization tables to be used in the encoding process
    /// @return The encoded image
    std::vector<unsigned char> encode_image(unsigned int Y, unsigned int X, const std::vector<Array_2d<double>>& arrays, 
        const std::vector<Comp_Info>& comp_infos, const std::vector<Huff_Table>& dc_tables, const std::vector<Huff_Table>& ac_tables,
        const std::vector<Q_Table>& q_tables);

    /// @brief Enlarges the component as necessary to ensure the required integer number of data units vertically and horizontally.
    /// If the original component satisfies these requriements, a copy is returned
    /// @param orig_comp Component to enlarge
    /// @param V Vertical sampling factor
    /// @param H Horizontal sampling factor
    /// @return A new Array_2d that satisfies the width and height requirements
    Array_2d<double> enlarge_component(const Array_2d<double>& orig_comp, unsigned int V, unsigned int H);

    /// @brief Transforms RGB to the YUV colour space in place as defined in ITU T.871
    /// @param red Red image component
    /// @param green Green image component
    /// @param blue Blue image component
    /// @return Tuple containing references to the YUV components
    std::tuple<Array_2d<double>&, Array_2d<double>&, Array_2d<double>&> colour_transform(
        Array_2d<double>& red, Array_2d<double>& green, Array_2d<double>& blue
    );
}

#endif
