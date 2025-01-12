#ifndef HUFF_TABLE_H
#define HUFF_TABLE_H

#include <cassert>
#include <vector>

#include "jpeg_cpp/general.h"
#include "jpeg_cpp/bit_string.h"

namespace JPEG
{
    class Huff_Table
    {
    private:
        std::vector<Bit_String> buffer;

    public:
        Huff_Table(size_t size) : buffer(size) {}

        Huff_Table(const std::initializer_list<Bit_String>& ls);

        size_t size() const { return buffer.size(); }

        Bit_String& operator[](size_t ind)
        {
            assert("Index out of bounds" && ind<size());

            return buffer[ind];
        }

        const Bit_String& operator[](size_t ind) const
        {
            assert("Index out of bounds" && ind<size());

            return buffer[ind];
        }

        static Huff_Table load_DC_table(Image_Component type);

        static Huff_Table load_AC_table(Image_Component type);
    };
}

#endif