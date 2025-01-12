#include "jpeg_cpp/component.h"

namespace JPEG
{
    Component::Component(Array_2d<double> array, size_t q_table_ind, size_t DC_Huff_table_ind, size_t AC_Huff_table_ind,
                            unsigned int H, unsigned int V)
    : array{ array }, q_table_ind{ q_table_ind }, DC_Huff_table_ind{ DC_Huff_table_ind },
        AC_Huff_table_ind{ AC_Huff_table_ind }, H{ H }, V{ V }
    {
        
    }
}

