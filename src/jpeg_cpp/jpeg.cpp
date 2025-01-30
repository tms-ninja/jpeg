#include "jpeg_cpp/jpeg.h"

namespace JPEG
{
    std::vector<unsigned char> encode_greyscale_image(const Array_2d<double>& array_2d)
    {
        // Load relevant quantization & Huffman tables
        std::vector<Q_Table> q_tables{
            Q_Table::load_spec_table(Image_Component::Luminance)
        };

        std::vector<Huff_Table> dc_huff_tables{
            Huff_Table::load_DC_table(Image_Component::Luminance)
        };

        std::vector<Huff_Table> ac_huff_tables{
            Huff_Table::load_AC_table(Image_Component::Luminance)
        };

        // Component info for this component
        std::vector<Comp_Info> comp_infos{
            Comp_Info{0, 0, 0, 1, 1}
        };

        // Enlarge component as neccessary
        std::vector<Array_2d<double>> enlarged_comps{
            enlarge_component(array_2d, comp_infos[0].V, comp_infos[0].H)
        };

        // Now call the general encode function
        const unsigned int Y{ static_cast<unsigned int>(array_2d.shape()[0]) };
        const unsigned int X{ static_cast<unsigned int>(array_2d.shape()[1]) };
        std::vector<unsigned char> encoded_image{ 
            encode_image(Y, X, enlarged_comps, comp_infos, dc_huff_tables, ac_huff_tables, q_tables) 
        };

        return encoded_image;
    }
}
