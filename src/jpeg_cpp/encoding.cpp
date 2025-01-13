#include "jpeg_cpp/encoding.h"

namespace JPEG
{
    void apply_level_shift(DU_Array<double> &array)
    {
        for (size_t ind = 0; ind < array.size(); ind++)
        {
            array[ind] -= 128.0;
        }
    }
}
