#include "jpeg_cpp/huff_table.h"

namespace JPEG
{
    Huff_Table::Huff_Table(const std::initializer_list<Bit_String> &ls)
    {
        for (auto& bs : ls)
        {
            buffer.push_back(bs);
        }
    }

    Huff_Table Huff_Table::load_DC_table(Image_Component type)
    {
        if (type == Image_Component::Luminance)
        {
            return Huff_Table{
                Bit_String{"00"},
                Bit_String{"010"},
                Bit_String{"011"},
                Bit_String{"100"},
                Bit_String{"101"},
                Bit_String{"110"},
                Bit_String{"1110"},
                Bit_String{"11110"},
                Bit_String{"111110"},
                Bit_String{"1111110"},
                Bit_String{"11111110"},
                Bit_String{"111111110"}
            };
        }

        // Return Chromiance table
        return Huff_Table{
            Bit_String{"00"},
            Bit_String{"01"},
            Bit_String{"10"},
            Bit_String{"110"},
            Bit_String{"1110"},
            Bit_String{"11110"},
            Bit_String{"111110"},
            Bit_String{"1111110"},
            Bit_String{"11111110"},
            Bit_String{"111111110"},
            Bit_String{"1111111110"},
            Bit_String{"11111111110"}
        };
    }

    Huff_Table Huff_Table::load_AC_table(Image_Component type)
    {
        Huff_Table table(0xFA+1);

        if (type==Image_Component::Luminance)
        {
            table[0x00] = Bit_String{"1010"};
            table[0x01] = Bit_String{"00"};
            table[0x02] = Bit_String{"01"};
            table[0x03] = Bit_String{"100"};
            table[0x04] = Bit_String{"1011"};
            table[0x05] = Bit_String{"11010"};
            table[0x06] = Bit_String{"1111000"};
            table[0x07] = Bit_String{"11111000"};
            table[0x08] = Bit_String{"1111110110"};
            table[0x09] = Bit_String{"1111111110000010"};
            table[0x0A] = Bit_String{"1111111110000011"};
            table[0x11] = Bit_String{"1100"};
            table[0x12] = Bit_String{"11011"};
            table[0x13] = Bit_String{"1111001"};
            table[0x14] = Bit_String{"111110110"};
            table[0x15] = Bit_String{"11111110110"};
            table[0x16] = Bit_String{"1111111110000100"};
            table[0x17] = Bit_String{"1111111110000101"};
            table[0x18] = Bit_String{"1111111110000110"};
            table[0x19] = Bit_String{"1111111110000111"};
            table[0x1A] = Bit_String{"1111111110001000"};
            table[0x21] = Bit_String{"11100"};
            table[0x22] = Bit_String{"11111001"};
            table[0x23] = Bit_String{"1111110111"};
            table[0x24] = Bit_String{"111111110100"};
            table[0x25] = Bit_String{"1111111110001001"};
            table[0x26] = Bit_String{"1111111110001010"};
            table[0x27] = Bit_String{"1111111110001011"};
            table[0x28] = Bit_String{"1111111110001100"};
            table[0x29] = Bit_String{"1111111110001101"};
            table[0x2A] = Bit_String{"1111111110001110"};
            table[0x31] = Bit_String{"111010"};
            table[0x32] = Bit_String{"111110111"};
            table[0x33] = Bit_String{"111111110101"};
            table[0x34] = Bit_String{"1111111110001111"};
            table[0x35] = Bit_String{"1111111110010000"};
            table[0x36] = Bit_String{"1111111110010001"};
            table[0x37] = Bit_String{"1111111110010010"};
            table[0x38] = Bit_String{"1111111110010011"};
            table[0x39] = Bit_String{"1111111110010100"};
            table[0x3A] = Bit_String{"1111111110010101"};
            table[0x41] = Bit_String{"111011"};
            table[0x42] = Bit_String{"1111111000"};
            table[0x43] = Bit_String{"1111111110010110"};
            table[0x44] = Bit_String{"1111111110010111"};
            table[0x45] = Bit_String{"1111111110011000"};
            table[0x46] = Bit_String{"1111111110011001"};
            table[0x47] = Bit_String{"1111111110011010"};
            table[0x48] = Bit_String{"1111111110011011"};
            table[0x49] = Bit_String{"1111111110011100"};
            table[0x4A] = Bit_String{"1111111110011101"};
            table[0x51] = Bit_String{"1111010"};
            table[0x52] = Bit_String{"11111110111"};
            table[0x53] = Bit_String{"1111111110011110"};
            table[0x54] = Bit_String{"1111111110011111"};
            table[0x55] = Bit_String{"1111111110100000"};
            table[0x56] = Bit_String{"1111111110100001"};
            table[0x57] = Bit_String{"1111111110100010"};
            table[0x58] = Bit_String{"1111111110100011"};
            table[0x59] = Bit_String{"1111111110100100"};
            table[0x5A] = Bit_String{"1111111110100101"};
            table[0x61] = Bit_String{"1111011"};
            table[0x62] = Bit_String{"111111110110"};
            table[0x63] = Bit_String{"1111111110100110"};
            table[0x64] = Bit_String{"1111111110100111"};
            table[0x65] = Bit_String{"1111111110101000"};
            table[0x66] = Bit_String{"1111111110101001"};
            table[0x67] = Bit_String{"1111111110101010"};
            table[0x68] = Bit_String{"1111111110101011"};
            table[0x69] = Bit_String{"1111111110101100"};
            table[0x6A] = Bit_String{"1111111110101101"};
            table[0x71] = Bit_String{"11111010"};
            table[0x72] = Bit_String{"111111110111"};
            table[0x73] = Bit_String{"1111111110101110"};
            table[0x74] = Bit_String{"1111111110101111"};
            table[0x75] = Bit_String{"1111111110110000"};
            table[0x76] = Bit_String{"1111111110110001"};
            table[0x77] = Bit_String{"1111111110110010"};
            table[0x78] = Bit_String{"1111111110110011"};
            table[0x79] = Bit_String{"1111111110110100"};
            table[0x7A] = Bit_String{"1111111110110101"};
            table[0x81] = Bit_String{"111111000"};
            table[0x82] = Bit_String{"111111111000000"};
            table[0x83] = Bit_String{"1111111110110110"};
            table[0x84] = Bit_String{"1111111110110111"};
            table[0x85] = Bit_String{"1111111110111000"};
            table[0x86] = Bit_String{"1111111110111001"};
            table[0x87] = Bit_String{"1111111110111010"};
            table[0x88] = Bit_String{"1111111110111011"};
            table[0x89] = Bit_String{"1111111110111100"};
            table[0x8A] = Bit_String{"1111111110111101"};
            table[0x91] = Bit_String{"111111001"};
            table[0x92] = Bit_String{"1111111110111110"};
            table[0x93] = Bit_String{"1111111110111111"};
            table[0x94] = Bit_String{"1111111111000000"};
            table[0x95] = Bit_String{"1111111111000001"};
            table[0x96] = Bit_String{"1111111111000010"};
            table[0x97] = Bit_String{"1111111111000011"};
            table[0x98] = Bit_String{"1111111111000100"};
            table[0x99] = Bit_String{"1111111111000101"};
            table[0x9A] = Bit_String{"1111111111000110"};
            table[0xA1] = Bit_String{"111111010"};
            table[0xA2] = Bit_String{"1111111111000111"};
            table[0xA3] = Bit_String{"1111111111001000"};
            table[0xA4] = Bit_String{"1111111111001001"};
            table[0xA5] = Bit_String{"1111111111001010"};
            table[0xA6] = Bit_String{"1111111111001011"};
            table[0xA7] = Bit_String{"1111111111001100"};
            table[0xA8] = Bit_String{"1111111111001101"};
            table[0xA9] = Bit_String{"1111111111001110"};
            table[0xAA] = Bit_String{"1111111111001111"};
            table[0xB1] = Bit_String{"1111111001"};
            table[0xB2] = Bit_String{"1111111111010000"};
            table[0xB3] = Bit_String{"1111111111010001"};
            table[0xB4] = Bit_String{"1111111111010010"};
            table[0xB5] = Bit_String{"1111111111010011"};
            table[0xB6] = Bit_String{"1111111111010100"};
            table[0xB7] = Bit_String{"1111111111010101"};
            table[0xB8] = Bit_String{"1111111111010110"};
            table[0xB9] = Bit_String{"1111111111010111"};
            table[0xBA] = Bit_String{"1111111111011000"};
            table[0xC1] = Bit_String{"1111111010"};
            table[0xC2] = Bit_String{"1111111111011001"};
            table[0xC3] = Bit_String{"1111111111011010"};
            table[0xC4] = Bit_String{"1111111111011011"};
            table[0xC5] = Bit_String{"1111111111011100"};
            table[0xC6] = Bit_String{"1111111111011101"};
            table[0xC7] = Bit_String{"1111111111011110"};
            table[0xC8] = Bit_String{"1111111111011111"};
            table[0xC9] = Bit_String{"1111111111100000"};
            table[0xCA] = Bit_String{"1111111111100001"};
            table[0xD1] = Bit_String{"11111111000"};
            table[0xD2] = Bit_String{"1111111111100010"};
            table[0xD3] = Bit_String{"1111111111100011"};
            table[0xD4] = Bit_String{"1111111111100100"};
            table[0xD5] = Bit_String{"1111111111100101"};
            table[0xD6] = Bit_String{"1111111111100110"};
            table[0xD7] = Bit_String{"1111111111100111"};
            table[0xD8] = Bit_String{"1111111111101000"};
            table[0xD9] = Bit_String{"1111111111101001"};
            table[0xDA] = Bit_String{"1111111111101010"};
            table[0xE1] = Bit_String{"1111111111101011"};
            table[0xE2] = Bit_String{"1111111111101100"};
            table[0xE3] = Bit_String{"1111111111101101"};
            table[0xE4] = Bit_String{"1111111111101110"};
            table[0xE5] = Bit_String{"1111111111101111"};
            table[0xE6] = Bit_String{"1111111111110000"};
            table[0xE7] = Bit_String{"1111111111110001"};
            table[0xE8] = Bit_String{"1111111111110010"};
            table[0xE9] = Bit_String{"1111111111110011"};
            table[0xEA] = Bit_String{"1111111111110100"};
            table[0xF0] = Bit_String{"11111111001"};
            table[0xF1] = Bit_String{"1111111111110101"};
            table[0xF2] = Bit_String{"1111111111110110"};
            table[0xF3] = Bit_String{"1111111111110111"};
            table[0xF4] = Bit_String{"1111111111111000"};
            table[0xF5] = Bit_String{"1111111111111001"};
            table[0xF6] = Bit_String{"1111111111111010"};
            table[0xF7] = Bit_String{"1111111111111011"};
            table[0xF8] = Bit_String{"1111111111111100"};
            table[0xF9] = Bit_String{"1111111111111101"};
            table[0xFA] = Bit_String{"1111111111111110"};

            return table;
        }

        table[0x00] = Bit_String{"00"};
        table[0x01] = Bit_String{"01"};
        table[0x02] = Bit_String{"100"};
        table[0x03] = Bit_String{"1010"};
        table[0x04] = Bit_String{"11000"};
        table[0x05] = Bit_String{"11001"};
        table[0x06] = Bit_String{"111000"};
        table[0x07] = Bit_String{"1111000"};
        table[0x08] = Bit_String{"111110100"};
        table[0x09] = Bit_String{"1111110110"};
        table[0x0A] = Bit_String{"111111110100"};
        table[0x11] = Bit_String{"1011"};
        table[0x12] = Bit_String{"111001"};
        table[0x13] = Bit_String{"11110110"};
        table[0x14] = Bit_String{"111110101"};
        table[0x15] = Bit_String{"11111110110"};
        table[0x16] = Bit_String{"111111110101"};
        table[0x17] = Bit_String{"1111111110001000"};
        table[0x18] = Bit_String{"1111111110001001"};
        table[0x19] = Bit_String{"1111111110001010"};
        table[0x1A] = Bit_String{"1111111110001011"};
        table[0x21] = Bit_String{"11010"};
        table[0x22] = Bit_String{"11110111"};
        table[0x23] = Bit_String{"1111110111"};
        table[0x24] = Bit_String{"111111110110"};
        table[0x25] = Bit_String{"111111111000010"};
        table[0x26] = Bit_String{"1111111110001100"};
        table[0x27] = Bit_String{"1111111110001101"};
        table[0x28] = Bit_String{"1111111110001110"};
        table[0x29] = Bit_String{"1111111110001111"};
        table[0x2A] = Bit_String{"1111111110010000"};
        table[0x31] = Bit_String{"11011"};
        table[0x32] = Bit_String{"11111000"};
        table[0x33] = Bit_String{"1111111000"};
        table[0x34] = Bit_String{"111111110111"};
        table[0x35] = Bit_String{"1111111110010001"};
        table[0x36] = Bit_String{"1111111110010010"};
        table[0x37] = Bit_String{"1111111110010011"};
        table[0x38] = Bit_String{"1111111110010100"};
        table[0x39] = Bit_String{"1111111110010101"};
        table[0x3A] = Bit_String{"1111111110010110"};
        table[0x41] = Bit_String{"111010"};
        table[0x42] = Bit_String{"111110110"};
        table[0x43] = Bit_String{"1111111110010111"};
        table[0x44] = Bit_String{"1111111110011000"};
        table[0x45] = Bit_String{"1111111110011001"};
        table[0x46] = Bit_String{"1111111110011010"};
        table[0x47] = Bit_String{"1111111110011011"};
        table[0x48] = Bit_String{"1111111110011100"};
        table[0x49] = Bit_String{"1111111110011101"};
        table[0x4A] = Bit_String{"1111111110011110"};
        table[0x51] = Bit_String{"111011"};
        table[0x52] = Bit_String{"1111111001"};
        table[0x53] = Bit_String{"1111111110011111"};
        table[0x54] = Bit_String{"1111111110100000"};
        table[0x55] = Bit_String{"1111111110100001"};
        table[0x56] = Bit_String{"1111111110100010"};
        table[0x57] = Bit_String{"1111111110100011"};
        table[0x58] = Bit_String{"1111111110100100"};
        table[0x59] = Bit_String{"1111111110100101"};
        table[0x5A] = Bit_String{"1111111110100110"};
        table[0x61] = Bit_String{"1111001"};
        table[0x62] = Bit_String{"11111110111"};
        table[0x63] = Bit_String{"1111111110100111"};
        table[0x64] = Bit_String{"1111111110101000"};
        table[0x65] = Bit_String{"1111111110101001"};
        table[0x66] = Bit_String{"1111111110101010"};
        table[0x67] = Bit_String{"1111111110101011"};
        table[0x68] = Bit_String{"1111111110101100"};
        table[0x69] = Bit_String{"1111111110101101"};
        table[0x6A] = Bit_String{"1111111110101110"};
        table[0x71] = Bit_String{"1111010"};
        table[0x72] = Bit_String{"11111111000"};
        table[0x73] = Bit_String{"1111111110101111"};
        table[0x74] = Bit_String{"1111111110110000"};
        table[0x75] = Bit_String{"1111111110110001"};
        table[0x76] = Bit_String{"1111111110110010"};
        table[0x77] = Bit_String{"1111111110110011"};
        table[0x78] = Bit_String{"1111111110110100"};
        table[0x79] = Bit_String{"1111111110110101"};
        table[0x7A] = Bit_String{"1111111110110110"};
        table[0x81] = Bit_String{"11111001"};
        table[0x82] = Bit_String{"1111111110110111"};
        table[0x83] = Bit_String{"1111111110111000"};
        table[0x84] = Bit_String{"1111111110111001"};
        table[0x85] = Bit_String{"1111111110111010"};
        table[0x86] = Bit_String{"1111111110111011"};
        table[0x87] = Bit_String{"1111111110111100"};
        table[0x88] = Bit_String{"1111111110111101"};
        table[0x89] = Bit_String{"1111111110111110"};
        table[0x8A] = Bit_String{"1111111110111111"};
        table[0x91] = Bit_String{"111110111"};
        table[0x92] = Bit_String{"1111111111000000"};
        table[0x93] = Bit_String{"1111111111000001"};
        table[0x94] = Bit_String{"1111111111000010"};
        table[0x95] = Bit_String{"1111111111000011"};
        table[0x96] = Bit_String{"1111111111000100"};
        table[0x97] = Bit_String{"1111111111000101"};
        table[0x98] = Bit_String{"1111111111000110"};
        table[0x99] = Bit_String{"1111111111000111"};
        table[0x9A] = Bit_String{"1111111111001000"};
        table[0xA1] = Bit_String{"111111000"};
        table[0xA2] = Bit_String{"1111111111001001"};
        table[0xA3] = Bit_String{"1111111111001010"};
        table[0xA4] = Bit_String{"1111111111001011"};
        table[0xA5] = Bit_String{"1111111111001100"};
        table[0xA6] = Bit_String{"1111111111001101"};
        table[0xA7] = Bit_String{"1111111111001110"};
        table[0xA8] = Bit_String{"1111111111001111"};
        table[0xA9] = Bit_String{"1111111111010000"};
        table[0xAA] = Bit_String{"1111111111010001"};
        table[0xB1] = Bit_String{"111111001"};
        table[0xB2] = Bit_String{"1111111111010010"};
        table[0xB3] = Bit_String{"1111111111010011"};
        table[0xB4] = Bit_String{"1111111111010100"};
        table[0xB5] = Bit_String{"1111111111010101"};
        table[0xB6] = Bit_String{"1111111111010110"};
        table[0xB7] = Bit_String{"1111111111010111"};
        table[0xB8] = Bit_String{"1111111111011000"};
        table[0xB9] = Bit_String{"1111111111011001"};
        table[0xBA] = Bit_String{"1111111111011010"};
        table[0xC1] = Bit_String{"111111010"};
        table[0xC2] = Bit_String{"1111111111011011"};
        table[0xC3] = Bit_String{"1111111111011100"};
        table[0xC4] = Bit_String{"1111111111011101"};
        table[0xC5] = Bit_String{"1111111111011110"};
        table[0xC6] = Bit_String{"1111111111011111"};
        table[0xC7] = Bit_String{"1111111111100000"};
        table[0xC8] = Bit_String{"1111111111100001"};
        table[0xC9] = Bit_String{"1111111111100010"};
        table[0xCA] = Bit_String{"1111111111100011"};
        table[0xD1] = Bit_String{"11111111001"};
        table[0xD2] = Bit_String{"1111111111100100"};
        table[0xD3] = Bit_String{"1111111111100101"};
        table[0xD4] = Bit_String{"1111111111100110"};
        table[0xD5] = Bit_String{"1111111111100111"};
        table[0xD6] = Bit_String{"1111111111101000"};
        table[0xD7] = Bit_String{"1111111111101001"};
        table[0xD8] = Bit_String{"1111111111101010"};
        table[0xD9] = Bit_String{"1111111111101011"};
        table[0xDA] = Bit_String{"1111111111101100"};
        table[0xE1] = Bit_String{"11111111100000"};
        table[0xE2] = Bit_String{"1111111111101101"};
        table[0xE3] = Bit_String{"1111111111101110"};
        table[0xE4] = Bit_String{"1111111111101111"};
        table[0xE5] = Bit_String{"1111111111110000"};
        table[0xE6] = Bit_String{"1111111111110001"};
        table[0xE7] = Bit_String{"1111111111110010"};
        table[0xE8] = Bit_String{"1111111111110011"};
        table[0xE9] = Bit_String{"1111111111110100"};
        table[0xEA] = Bit_String{"1111111111110101"};
        table[0xF0] = Bit_String{"1111111010"};
        table[0xF1] = Bit_String{"111111111000011"};
        table[0xF2] = Bit_String{"1111111111110110"};
        table[0xF3] = Bit_String{"1111111111110111"};
        table[0xF4] = Bit_String{"1111111111111000"};
        table[0xF5] = Bit_String{"1111111111111001"};
        table[0xF6] = Bit_String{"1111111111111010"};
        table[0xF7] = Bit_String{"1111111111111011"};
        table[0xF8] = Bit_String{"1111111111111100"};
        table[0xF9] = Bit_String{"1111111111111101"};
        table[0xFA] = Bit_String{"1111111111111110"};

        return table;
    }
}
