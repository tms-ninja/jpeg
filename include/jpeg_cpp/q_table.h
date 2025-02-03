#ifndef Q_TABLE_H
#define Q_TABLE_H

#include <array>
#include <cmath>

#include "jpeg_cpp/array.h"
#include "general.h"

namespace JPEG
{
    /// @brief Represents a quantization table
    class Q_Table
    {
    private:
        Array_2d<unsigned short> arr_2d;

    public:
        /// @brief Constructs an uninitialized Q_Table
        Q_Table() : arr_2d{8, 8} {}

        /// @brief Constructs a Q_Table from nested initialization lists
        /// @param ls Initializer list with shape (8, 8)
        Q_Table(const std::initializer_list<std::initializer_list<unsigned short>>& ls);

        /// @brief Pointer to the buffer underlying the Q_Table
        /// @return Pointer to the buffer underlying the Q_Table
        const unsigned short* data() const { return arr_2d.data(); }

        /// @brief Total number of elements in the Q_Table
        /// @return Total number of elements in the Q_Table
        size_t size() const { return arr_2d.size(); }

        /// @brief Shape of the Q_Table
        /// @return Shape of the Q_Table
        std::array<size_t, 2> shape() const { return arr_2d.shape(); }

        /// @brief Access elements as a 1d array
        /// @param ind Index of element as a 1d array
        /// @return Element at index ind
        unsigned short& operator[](size_t ind) { return arr_2d[ind]; }

        /// @brief Access elements as a 1d array
        /// @param ind Index of element as a 1d array
        /// @return Element at index ind
        const unsigned short& operator[](size_t ind) const { return arr_2d[ind]; }

        /// @brief Gets the element at the given row and column
        /// @param row The row of the desired element
        /// @param column The column of the desired element
        /// @return Reference to the element at the given position
        unsigned short& at(size_t row, size_t column) { return arr_2d.at(row, column); }
        
        /// @brief Gets the element at the given row and column
        /// @param row The row of the desired element
        /// @param column The column of the desired element
        /// @return Reference to the element at the given position
        const unsigned short& at(size_t row, size_t column) const { return arr_2d.at(row, column); }

        /// @brief Gets the element at the given row and column
        /// @param row The row of the desired element
        /// @param column The column of the desired element
        /// @return Reference to the element at the given position
        unsigned short& operator()(size_t row, size_t column) { return arr_2d.at(row, column); }

        /// @brief Gets the element at the given row and column
        /// @param row The row of the desired element
        /// @param column The column of the desired element
        /// @return Reference to the element at the given position
        const unsigned short& operator()(size_t row, size_t column) const { return arr_2d.at(row, column); }

        friend std::ostream& operator<<(std::ostream& out, const Q_Table& q_table);

        /// @brief Loads an example quantization table from the JPEG spec
        /// @param type Component type (Luminance/Chromiance)
        /// @return Quantization table
        static Q_Table load_spec_table(Image_Component type);

        /// @brief Generates a quantization table for a given quality factor using the IJG algorithm.
        /// @param type Component type (Luminance/Chromiance)
        /// @param qf Quality factor, integer from 0 to 100 inclusive. 50 corresponds to JPEG spec suggested
        /// tables
        /// @return Quantization table
        static Q_Table load_q_table_from_quality_factor(Image_Component type, int qf);
    };
}

#endif
