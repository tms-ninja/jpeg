#ifndef ARRAY_H
#define ARRAY_H

#include <array>
#include <iostream>
#include <vector>

namespace JPEG {

    /// @brief Represents a contiguous array of elements of type T
    /// @tparam T Element type
    template<typename T>
    class Array
    {
    private:
        std::vector<T> buffer;

    public:
        /// @brief Creates an empty Array
        Array() {}

        /// @brief Creates an Array of the given size
        /// @param size Number of elements to create array with
        Array(size_t size) : buffer(size) {}

        /// @brief Initializes Array with given elements
        /// @param ls Initializer list containing the elements 
        Array(const std::initializer_list<T>& ls)
        {
            for (auto& elem: ls)
            {
                buffer.push_back(elem);
            }
        }

        /// @brief Number of elements in the Array
        /// @return Number of elements in the Array
        size_t size() const { return buffer.size(); }

        /// @brief Returns a pointer to the underlying buffer
        /// @return A pointer to the underlying buffer
        const T* data() const { return buffer.data(); }

        /// @brief Resizes the array
        /// @param size New size of the array
        void resize(size_t size)
        {
            buffer.resize(size);
        }

        T& operator[](size_t ind) { return buffer[ind] ;}
        const T& operator[](size_t ind) const { return buffer[ind] ;}

        friend std::ostream& operator<< (std::ostream& out, const Array<T>& arr)
        {
            for (T& elem: arr.buffer)
            {
                out << elem << ' ';
            }

            return out;
        }
    };

    /// @brief A 2d array stored as a contiguous array in row-major order
    /// @tparam T Element type
    template<typename T>
    class Array_2d
    {
    private:
        Array<T> buffer;
        std::array<size_t, 2> arr_shape;

    public:
        /// @brief Creates an Array_2d with the given number of rows and columns
        /// @param rows 
        /// @param columns 
        Array_2d(size_t rows, size_t columns)
        : buffer(rows*columns), arr_shape{rows, columns}
        {

        }

        Array_2d(const std::initializer_list<std::initializer_list<T>>& ls)
        {
            // For now, assume neither number of rows or columns is zero
            arr_shape[0] = ls.size();
            arr_shape[1] = (*ls.begin()).size();

            buffer.resize(arr_shape[0]*arr_shape[1]);

            size_t i{ 0 }, j{ 0 };

            for (auto& row: ls)
            {
                // Verify the row size

                // remember to reset j
                j = 0;
                
                for (auto& elem: row)
                {
                    at(i, j) = elem;
                    j++;
                }

                i++;
            }
        }

        /// @brief Shape of the array
        /// @return Shape of the array as rows/columns
        std::array<size_t, 2> shape() const { return arr_shape; }

        /// @brief Number of elements in the Array_2d
        /// @return Number of elements in the Array_2d
        size_t size() const { return buffer.size(); }

        /// @brief Returns a pointer to the underlying buffer
        /// @return A pointer to the underlying buffer
        const T* data() const { return buffer.data(); }

        /// @brief Access elements as a 1d array
        /// @param ind Index of element as a 1d array
        /// @return Element at index ind
        T& operator[](size_t ind) { return buffer[ind] ;}

        /// @brief Access elements as a 1d array
        /// @param ind Index of element as a 1d array
        /// @return Element at index ind
        const T& operator[](size_t ind) const { return buffer[ind] ;}

        /// @brief Gets the element at the given row and column
        /// @param row The row of the desired element
        /// @param column The column of the desired element
        /// @return Reference to the element at the given position
        T& at(size_t row, size_t column)
        {
            size_t ind{ row*arr_shape[1] + column };

            return buffer[ind];
        }
        
        /// @brief Gets the element at the given row and column
        /// @param row The row of the desired element
        /// @param column The column of the desired element
        /// @return Reference to the element at the given position
        const T& at(size_t row, size_t column) const
        {
            size_t ind{ row*arr_shape[1] + column };

            return buffer[ind];
        }

        /// @brief Gets the element at the given row and column
        /// @param row The row of the desired element
        /// @param column The column of the desired element
        /// @return Reference to the element at the given position
        T& operator()(size_t row, size_t column)
        {
            return at(row, column);
        }

        /// @brief Gets the element at the given row and column
        /// @param row The row of the desired element
        /// @param column The column of the desired element
        /// @return Reference to the element at the given position
        const T& operator()(size_t row, size_t column) const
        {
            return at(row, column);
        }

        friend std::ostream& operator<<(std::ostream& out, const Array_2d<T>& arr)
        {
            for (size_t i=0; i<arr.arr_shape[0]; i++)
            {
                for (size_t j=0; j<arr.arr_shape[1]-1; j++)
                {
                    out << at(i, j) << ' ';
                }
                // Don't forget the last element in the row
                out << at(i, arr.arr_shape[1]-1) << '\n';
            }
            
            return out;
        }
    };
}

#endif
