#include "jpeg_cpp/q_table.h"

namespace JPEG
{
    Q_Table::Q_Table(const std::initializer_list<std::initializer_list<unsigned short>> &ls)
     : arr_2d{ 8, 8 }
    {
        // Need to check the ls is 8 by 8
        assert("Q_Tables must have shape (8, 8)" && ls.size()==8);

        for (auto& row: ls)
        {
            assert("Q_Tables must have shape (8, 8)" && row.size()==8);
        }

        // All good so copy to the array
        arr_2d = Array_2d<unsigned short>{ls};
    }

    std::ostream& operator<<(std::ostream& out, const Q_Table& q_table)
    {
        out << q_table.arr_2d;

        return out;
    }

    Q_Table Q_Table::load_spec_table(Image_Component type)
    {
        Q_Table q_table;

        if (type==Image_Component::Luminance)
        {
            q_table = Q_Table{
                {16, 11, 10, 16,  24,  40,  51,  61},
                {12, 12, 14, 19,  26,  58,  60,  55},
                {14, 13, 16, 24,  40,  57,  69,  56},
                {14, 17, 22, 29,  51,  87,  80,  62},
                {18, 22, 37, 56,  68, 109, 103,  77},
                {24, 35, 55, 64,  81, 104, 113,  92},
                {49, 64, 78, 87, 103, 121, 120, 101},
                {72, 92, 95, 98, 112, 100, 103,  99},
            };
        }
        else
        {
            q_table = Q_Table{
                {17, 18, 24, 47, 99, 99, 99, 99},
                {18, 21, 26, 66, 99, 99, 99, 99},
                {24, 26, 56, 99, 99, 99, 99, 99},
                {47, 66, 99, 99, 99, 99, 99, 99},
                {99, 99, 99, 99, 99, 99, 99, 99},
                {99, 99, 99, 99, 99, 99, 99, 99},
                {99, 99, 99, 99, 99, 99, 99, 99},
                {99, 99, 99, 99, 99, 99, 99, 99},
            };
        }

        return q_table;
    }

    Q_Table Q_Table::load_q_table_from_quality_factor(Image_Component type, int qf)
    {
        if (qf<0)
        {
            throw std::invalid_argument("Quality factor cannot be less than 0");
        }
        else if (qf>100)
        {
            throw std::invalid_argument("Quality factor cannot be greater than 100");
        }

        Q_Table table{ Q_Table::load_spec_table(type) };
        
        // Avoid division of zero later
        if (qf==0)
        {
            qf = 1;
        }

        double s{ qf<50 ? 5000.0/qf : 200.0-2.0*qf };

        for (size_t ind = 0; ind < table.size(); ind++)
        {
            table[ind] = std::floor(s*table[ind]/100.0 + 0.5);

            if (table[ind]==0)
            {
                table[ind] = 1;
            }
            else if (table[ind]>255)
            {
                // quantization table elements of 8 bit precision images
                // can't be more than 8 bits
                table[ind] = 255;
            }
        }

        return table;
    }
}