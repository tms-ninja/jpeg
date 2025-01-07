#ifndef ARRAY_H
#define ARRAY_H

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
}

#endif
