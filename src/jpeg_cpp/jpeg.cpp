#include "jpeg_cpp/jpeg.h"

namespace JPEG
{
    std::vector<unsigned char> encode_greyscale_image(const Array_2d<double>& array_2d, int qf)
    {
        // Load relevant quantization & Huffman tables
        std::vector<Q_Table> q_tables{
            Q_Table::load_q_table_from_quality_factor(Image_Component::Luminance, qf)
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
        const Array_2d<double>& red, const Array_2d<double>& green, const Array_2d<double>& blue, int qf,
        Subsampling ss
    )
    {
        // Load relevant quantization & Huffman tables
        std::vector<Q_Table> q_tables{
            Q_Table::load_q_table_from_quality_factor(Image_Component::Luminance, qf),
            Q_Table::load_q_table_from_quality_factor(Image_Component::Chrominance, qf)
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
        std::vector<Comp_Info> comp_infos;

        switch (ss)
        {
        case Subsampling::ss_4_4_4:
            // No subsampling, H and V both 1
            comp_infos.emplace_back(0, 0, 0, 1, 1);
            break;
        case Subsampling::ss_4_2_2:
            comp_infos.emplace_back(0, 0, 0, 2, 1);
            break;
        case Subsampling::ss_4_2_0:
            comp_infos.emplace_back(0, 0, 0, 2, 2);
            break;
        default:
            throw std::invalid_argument("Chosen subsampling is not supported");
        }

        // Chromiance components aways have horizontal and vertical subsampling factors set to 1
        comp_infos.emplace_back(1, 1, 1, 1, 1);
        comp_infos.emplace_back(1, 1, 1, 1, 1);

        // Enlarge component as neccessary. Since H and V are always either 1 or 2 (in our case)
        // we just need to find the maximum H and V for all components
        unsigned int max_H{
            std::max_element(comp_infos.begin(), comp_infos.end(), 
                [](Comp_Info a, Comp_Info b)
                {
                    return a.H < b.H;
                }
            )->H
        };
        unsigned int max_V{
            std::max_element(comp_infos.begin(), comp_infos.end(), 
                [](Comp_Info a, Comp_Info b)
                {
                    return a.V < b.V;
                }
            )->V
        };

        std::vector<Array_2d<double>> enlarged_comps{
            enlarge_component(red, max_V, max_H),
            enlarge_component(green, max_V, max_H),
            enlarge_component(blue, max_V, max_H)
        };

        // Perform the colour transform
        colour_transform(enlarged_comps[0], enlarged_comps[1], enlarged_comps[2]);

        // Only subsample the chromiance components
        for (size_t ind = 1; ind < enlarged_comps.size(); ind++)
        {
            auto &comp{ enlarged_comps[ind] };
            subsample_component(comp, ss);
        }

        // Now call the general encode function
        const unsigned int Y{ static_cast<unsigned int>(red.shape()[0]) };
        const unsigned int X{ static_cast<unsigned int>(red.shape()[1]) };
        std::vector<unsigned char> encoded_image{ 
            encode_image(Y, X, enlarged_comps, comp_infos, dc_huff_tables, ac_huff_tables, q_tables) 
        };

        return encoded_image;
    }
}
