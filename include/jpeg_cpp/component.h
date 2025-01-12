#ifndef COMPONENT_H
#define COMPONENT_H

#include "jpeg_cpp/array.h"

namespace JPEG
{   
    /// @brief Represents an image component.
    struct Component
    {
        /// @brief Image component data
        Array_2d<double> array;

        /// @brief Index of the Quantization table for this component
        size_t q_table_ind;

        /// @brief Index of the DC Huffman table for this component
        size_t DC_Huff_table_ind;

        /// @brief Index of the AC Huffman table for this component
        size_t AC_Huff_table_ind;
        
        /// @brief Horizontal sampling factor
        unsigned int H;

        /// @brief Vertical sampling factor
        unsigned int V;

        /// @brief Constructs a Component
        /// @param array Array containing component data
        /// @param q_table_ind Index of the Quantization table
        /// @param DC_Huff_table_ind Index of the DC Huffman table
        /// @param AC_Huff_table_ind Index of the AC Huffman table
        /// @param H Horizontal sampling factor
        /// @param V Vertical sampling factor
        Component(Array_2d<double> array, size_t q_table_ind, size_t DC_Huff_table_ind, size_t AC_Huff_table_ind,
                    unsigned int H, unsigned int V);
    };
    
}

#endif
