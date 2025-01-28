#ifndef HUFF_TABLE_H
#define HUFF_TABLE_H

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

        /// @brief 
        /// @param type 
        /// @return 
        static Huff_Table load_AC_table(Image_Component type);
    };
}

#endif