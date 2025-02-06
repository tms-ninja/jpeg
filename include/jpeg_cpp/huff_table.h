#ifndef HUFF_TABLE_H
#define HUFF_TABLE_H

#include <algorithm>
#include <array>
#include <cassert>
#include <vector>

#include "jpeg_cpp/general.h"
#include "jpeg_cpp/bit_string.h"

namespace JPEG
{
    /// @brief Represents a Huffman table
    class Huff_Table
    {
    private:
        std::vector<Bit_String> buffer;

        /// @brief Determines the code size for Huffman codes using the frequency data. Follows the
        /// algorithm in Figure K.1 of the JPEG spec
        /// @param freq Frequency data of symbols
        /// @param code_size Code size for Huffman codes, should be initialized to zero
        /// @param others The others array
        static void compute_code_size(std::array<unsigned int, 257>& freq, std::array<unsigned int, 257>& code_size, std::array<int, 257>& others);

        /// @brief Computes the BITS array given how many codes there are for a given code length.
        /// Follows the algorithm of Figure K.2 of the JPEG spec
        /// @param code_size Code size of each symbol
        /// @return BITS array
        static std::array<unsigned int, 16> count_BITS(const std::array<unsigned int, 257>& code_size);

        /// @brief Adjusts the code ssizes so that no code has more than 16 bits. Follows the algorithm
        /// in Figure K.3 of the JPEG spec
        /// @param bits the BITS array
        static void adjust_BITS(std::array<unsigned int, 32>& bits);
        
        /// @brief Generates the HUFFVAL array, follows the algorithm in Figure K.4 of the JPEG spec
        /// @param code_size_array Vector of code sizes for each symbol
        /// @return HUFFVAL array
        static std::vector<unsigned int> sort_input(const std::array<unsigned int, 257>& code_size_array);
    public:
        /// @brief Constructs a Huff_Table with the given number of elements. The Huffman
        /// codes of the newly constructed elements initially have size 0.
        /// @param size Number of entries in the Huff_Table
        Huff_Table(size_t size) : buffer(size) {}

        /// @brief Constructs a Huff_Table from the given Huffman codes.
        /// @param ls List of Huffman codes
        Huff_Table(const std::initializer_list<Bit_String>& ls);

        /// @brief Number of entries in the Huffman table
        /// @return Total number of entries in the table including codes of zero size
        size_t size() const { return buffer.size(); }

        /// @brief Gets the Huffman code corresponding to the given index
        /// @param ind Desired index of the Huffman code
        /// @return Huffman code corresponding to the given index
        Bit_String& operator[](size_t ind)
        {
            assert("Index out of bounds" && ind<size());

            return buffer[ind];
        }

        /// @brief Gets the Huffman code corresponding to the given index
        /// @param ind Desired index of the Huffman code
        /// @return Huffman code corresponding to the given index
        const Bit_String& operator[](size_t ind) const
        {
            assert("Index out of bounds" && ind<size());

            return buffer[ind];
        }

        /// @brief Returns a spec DC Huffman table
        /// @param type Type of component the table corresponds to, luminance or chromiance
        /// @return Spec DC Huffman table
        static Huff_Table load_DC_table(Image_Component type);

        /// @brief Returns a spec AC Huffman table
        /// @param type Type of component the table corresponds to, luminance or chromiance
        /// @return Spec AC Huffman table
        static Huff_Table load_AC_table(Image_Component type);

        /// @brief Loads the Huffman table from the given BITS and HUFFVAL arrays
        /// @param bits BITS array
        /// @param huffval HUFFVAL array
        /// @return Huffman table loaded from the given BITS and HUFFVAL arrays
        static Huff_Table load_table_from_BITS_and_HUFFVAL(
            const std::array<unsigned int, 16>& bits, const std::vector<unsigned int>& huffval
        );

        /// @brief Generates an optimal Huffman table for the given stats. Number of symbols 
        /// should not exceed 256
        /// @param stats Frequency count of each symbol
        /// @return Optimal Huffman table for the given frequencies
        static Huff_Table gen_table_from_stats(const std::vector<unsigned int>& stats);
    };
}

#endif