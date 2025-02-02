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

    std::vector<unsigned char> encode_colour_image(
        const Array_2d<double>& red, const Array_2d<double>& green, const Array_2d<double>& blue
    )
    {
        // Load relevant quantization & Huffman tables
        std::vector<Q_Table> q_tables{
            Q_Table::load_spec_table(Image_Component::Luminance),
            Q_Table::load_spec_table(Image_Component::Chrominance)
        };

        std::vector<Huff_Table> dc_huff_tables{
            Huff_Table::load_DC_table(Image_Component::Luminance),
            Huff_Table::load_DC_table(Image_Component::Chrominance)
        };

        std::vector<Huff_Table> ac_huff_tables{
            Huff_Table::load_AC_table(Image_Component::Luminance),
            Huff_Table::load_AC_table(Image_Component::Chrominance)
        };

        // Component info for the components
        // Arguments are:
        // q_table_ind
        // DC_Huff_table_ind
        // AC_Huff_table_ind
        // H
        // V
        std::vector<Comp_Info> comp_infos{
            // Luminance component
            Comp_Info{0, 0, 0, 1, 1},
            // Chromiance components
            // We'll ignore subsampling for now
            Comp_Info{1, 1, 1, 1, 1},
            Comp_Info{1, 1, 1, 1, 1},
        };

        // Enlarge component as neccessary
        std::vector<Array_2d<double>> enlarged_comps{
            // enlarge_component() needs to be updated to use number of samples 
            // instead of sampling factors  
            enlarge_component(red, comp_infos[0].V, comp_infos[0].H),
            enlarge_component(green, comp_infos[1].V, comp_infos[1].H),
            enlarge_component(blue, comp_infos[2].V, comp_infos[2].H)
        };

        // Perform the colour transform
        colour_transform(enlarged_comps[0], enlarged_comps[1], enlarged_comps[2]);

        // If future we'll perform subsampling here but we'll leave that for now

        // Now call the general encode function
        const unsigned int Y{ static_cast<unsigned int>(red.shape()[0]) };
        const unsigned int X{ static_cast<unsigned int>(red.shape()[1]) };
        std::vector<unsigned char> encoded_image{ 
            encode_image(Y, X, enlarged_comps, comp_infos, dc_huff_tables, ac_huff_tables, q_tables) 
        };

        return encoded_image;
    }
}
