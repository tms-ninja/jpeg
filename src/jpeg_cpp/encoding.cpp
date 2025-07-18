#include "jpeg_cpp/encoding.h"

namespace JPEG
{
    void append_two_bytes(std::vector<unsigned char>& out, unsigned long long bytes)
    {
        out.push_back(static_cast<unsigned char>(bytes >> 8));
        out.push_back(static_cast<unsigned char>(0xFF & bytes));
    }

    void append_composite_byte(std::vector<unsigned char>& out, unsigned long long high, unsigned long long low)
    {
        out.push_back(static_cast<unsigned char>((high << 4) | low));
    }

    void append_marker(std::vector<unsigned char>& out, Marker marker)
    {
        append_two_bytes(out, static_cast<unsigned long long>(marker));
    }

    void apply_level_shift(DU_Array<double> &array)
    {
        for (size_t ind = 0; ind < array.size(); ind++)
        {
            array[ind] -= 128.0;
        }
    }

    Array_2d<double> load_DCT_matrix()
    {
        Array_2d<double> dct_matrix{ 8, 8 };

        for (size_t i = 0; i < dct_matrix.shape()[0]; i++)
        {
            for (size_t j = 0; j < dct_matrix.shape()[1]; j++)
            {
                dct_matrix(i, j) = std::cos((2*i+1)*j*M_PI / 16.0) / 2.0;
            }
        }
        
        return dct_matrix;
    }

    Array_2d<double> load_DCT_matrix_transpose()
    {
        Array_2d<double> dct_matrix{ load_DCT_matrix() };

        // Compute the transpose manually
        for (size_t i = 0; i < dct_matrix.shape()[0]; i++)
        {
            for (size_t j = 0; j < i; j++)
            {
                double old_value_at_ij{ dct_matrix(i, j) };

                dct_matrix(i, j) = dct_matrix(j, i);
                dct_matrix(j, i) = old_value_at_ij;
            }
        }

        return dct_matrix;
    }

    void mat_mul(const double* a, const double* b, double* c, size_t N)
    {
        for (size_t i = 0; i < N; i++)
        {
            for (size_t j = 0; j < N; j++)
            {
                // Compute c_ij = sum_k a_ik*b_kj
                double c_ij{ 0.0 };

                for (size_t k = 0; k < N; k++)
                {
                    c_ij += a[i*N+k] * b[k*N+j];
                }
                
                // Don't forget to set c_ij!
                c[i*N + j] = c_ij;
            }
        }
    }

    void apply_DCT(DU_Array<double> &array)
    {
        // Buffer to store the intermediate result of data unit*DCT matrix
        std::array<double, du_size> double_buffer;
        const auto arr_shape{ array.shape() };

        for (size_t du_ind = 0; du_ind < arr_shape[0]; du_ind++)
        {
            // Compute data unit*DCT matrix
            double* du_ptr{ &array[du_ind*du_size] };

            mat_mul(du_ptr, DCT_matrix.data(), &double_buffer[0], arr_shape[1]);

            // Compute transpose(DCT matrix) * termporary
            mat_mul(DCT_matrix_transpose.data(), &double_buffer[0], du_ptr, arr_shape[1]);
        }
    }

    void apply_quantization(DU_Array<double> &array, const Q_Table& q_table)
    {
        // Need to apply factors of sqrt(2) to the elements of the q_table that
        // come from the DCT
        Array_2d<double> q_table_corrected{ q_table.shape()[0],  q_table.shape()[1] };

        for (size_t i = 0; i < q_table_corrected.shape()[0]; i++)
        {
            for (size_t j = 0; j < q_table_corrected.shape()[1]; j++)
            {
                q_table_corrected(i, j) = 1.0/q_table(i, j);

                if (i==0)
                {
                    q_table_corrected(i, j) /= std::sqrt(2.0);
                }
                
                if (j==0)
                {
                    q_table_corrected(i, j) /= std::sqrt(2.0);
                }
            }
        }

        size_t data_unit_offset{ 0 };

        // Now we can apply quantization
        for (size_t du_ind = 0; du_ind < array.shape()[0]; du_ind++)
        {

            for (size_t ind = 0; ind < du_size; ind++)
            {
                array[data_unit_offset+ind] = std::round(array[data_unit_offset+ind] * q_table_corrected[ind]);
            }

            data_unit_offset += du_size;
        }
    }

    unsigned int compute_ssss(unsigned int n)
    {
        // Check n<=2047 as expected
        if (n>2047)
        {
            throw std::invalid_argument("n cannot be greater than 2047");
        }

        if (n==0)
        {
            return 0;
        }

        // Assume n<=2047 so that ssss<=11
        unsigned int ssss{ 11 };

        while ((n >> ssss) == 0)
        {
            ssss -= 1;
        }

        return ssss+1;
    }

    void encode_DC_coeff(std::vector<Coefficient>& coeffs, int diff, size_t comp_ind)
    {
        unsigned int ssss{ compute_ssss(std::abs(diff)) };

        coeffs.emplace_back(diff, ssss, Coefficient_Type::DC, comp_ind);
    }
    
    void encode_AC_coeff(std::vector<Coefficient>& coeffs, int ac, unsigned int rrrr, size_t comp_ind)
    {
        assert("coeff should not be zero" && ac!=0);
        assert("rrrr must be less than 16" && rrrr<16);

        unsigned int ssss{ compute_ssss(std::abs(ac)) };
        unsigned int rrrrssss{ (rrrr << 4) + ssss };

        coeffs.emplace_back(ac, rrrrssss, Coefficient_Type::AC, comp_ind);
    }

    void encode_AC_coeffs(std::vector<Coefficient>& coeffs, const DU_Array<double>& du_array, size_t du_ind, size_t comp_ind)
    {
        // Offsets for the zig-zag pattern
        const auto& zz{ zig_zag_order };

        // Pointer to the start of the data unit we are encoding
        const double* du_ptr{ du_array.data() + du_ind*du_size};
        size_t end_of_block{ 0x00 };  // EOB RRRRSSSS code
        size_t zrl{ 0xF0 };  // ZRL RRRRSSSS code (run of 16 zeros)

        // Index of the current coefficient we are encoding. Not k=0 correpsonds
        // to the DC coefficient that we do not encode in this function
        size_t k{ 0 };
        size_t r{ 0 };  // Current run length of zeros. Can be much higher than 16, i.e. 32+

        while (k<63)
        {
            k++;

            // Get the next value to encode
            int value_to_encode{ static_cast<int>(*(du_ptr+zz[k])) };

            if (value_to_encode==0)
            {
                if (k==63)
                {
                    // No more values to encode, encode EOB
                    // Use dummy value of zero as note bits from it are encoded
                    coeffs.emplace_back(0, end_of_block, Coefficient_Type::AC, comp_ind);
                }

                r++;
                continue;
            }

            // Note r could be much larger than 15 and might need to encode more 
            // than one ZRL
            while (r>15)
            {
                // Encode ZRL
                // Use dummy value of zero as note bits from it are encoded
                coeffs.emplace_back(0, zrl, Coefficient_Type::AC, comp_ind);
                r -= 16;
            }

            // Now encode the current, non-zero coefficient
            encode_AC_coeff(coeffs, value_to_encode, static_cast<unsigned char>(r), comp_ind);

            // Reset the run length of zeros
            r = 0;
        }
    }

    int encode_data_unit_sequential(std::vector<Coefficient>& coeffs, const DU_Array<double>& du_array, size_t du_ind, int prev_dc,
        size_t comp_ind)
    {
        // First encode the DC coefficient
        int dc_coeff{ static_cast<int>(du_array(du_ind, 0, 0)) };
        // Its the difference between the DC coeff in this unit and the previous one that
        // is encoded
        int dc_diff{ dc_coeff-prev_dc };

        // encode_DC_coeff(bs, dc_diff, huff_table_dc);
        encode_DC_coeff(coeffs, dc_diff, comp_ind);

        // Now encode the AC coefficients
        // encode_AC_coeffs(bs, du_array, du_ind, huff_table_ac);
        encode_AC_coeffs(coeffs, du_array, du_ind, comp_ind);

        return dc_coeff;
    }

    void append_q_table_marker_segment(std::vector<unsigned char>& out, const std::vector<Q_Table>& q_tables)
    {
        // Append quantization table marker, 0xFFDB, most significant byte first
        append_marker(out, Marker::Define_Quantization_Table);

        // Now determine the length of the marker segment
        size_t length{ 2 + 65 * q_tables.size() };

        append_two_bytes(out, length);

        // Now append each table in turn
        for (size_t table_ind = 0; table_ind < q_tables.size(); table_ind++)
        {
            const Q_Table& table{ q_tables[table_ind] };

            // Table precision, 0 for 8 bit images
            unsigned int Pq{ 0 };

            // Table destination identifier
            unsigned int Tq{ static_cast<unsigned int>(table_ind) };
            assert("Tq cannot be greater than 3" && Tq<=3);

            // Pq & Tq form one byte, Pq most significant
            append_composite_byte(out, Pq, Tq);

            // Now need to append the value of the quantization table in zig-zag order
            const auto& zz{ zig_zag_order };

            for (size_t ind = 0; ind < table.size(); ind++)
            {   
                assert("Quantization table values must be less than 256" && table[zz[ind]]);
                out.push_back(static_cast<unsigned char>(table[zz[ind]]));
            }
        }
    }

    void append_mcu(std::vector<Coefficient>& coeffs, std::vector<int>& prev_dc, std::vector<size_t>& du_ind, const std::vector<DU_Array<double>>& arrays, 
        const std::vector<Comp_Info>& comp_infos)
    {
        // Iterate over each component adding however many data units we need from it
        for (size_t comp_ind = 0; comp_ind < arrays.size(); comp_ind++)
        {
            const DU_Array<double>& data_units{ arrays[comp_ind] };

            size_t data_units_to_encode{ comp_infos[comp_ind].H*comp_infos[comp_ind].V };

            for (size_t ind = 0; ind < data_units_to_encode; ind++)
            {
                // encode the current data unit and update the current
                // DC coeff for the current component
                prev_dc[comp_ind] = encode_data_unit_sequential(
                    coeffs, 
                    data_units, 
                    du_ind[comp_ind],
                    prev_dc[comp_ind],
                    comp_ind
                );

                // Move on to next data unit
                du_ind[comp_ind]++;
            }
        }
    }

    std::vector<Coefficient> encode_coeff_rep_sequential(const std::vector<DU_Array<double>>& arrays, const std::vector<Comp_Info>& comp_infos)
    {
        // Entropy encode each MCU in turn
        // Note previous DC should be initialized to 0
        std::vector<int> prev_dc(arrays.size());
        std::vector<size_t> du_ind(arrays.size());

        std::vector<Coefficient> encoded_coeffs;

        while (du_ind[0]<arrays[0].shape()[0])
        {
            append_mcu(encoded_coeffs, prev_dc, du_ind, arrays, comp_infos);         
        }

        return encoded_coeffs;
    }

    void append_BITS_array(std::vector<unsigned char>& out, const Huff_Table& huff_table)
    {
        std::array<unsigned char, 16> BITS{};

        for (size_t ind = 0; ind < huff_table.size(); ind++)
        {
            const Bit_String& bs{ huff_table[ind] };

            assert("Bit_String length must be less than 17" && bs.size()<17);

            if (bs.size()>0)
            {
                BITS[bs.size()-1]++;
            }
        }
        
        for (auto b : BITS)
        {
            out.push_back(b);
        }
    }

    void append_HUFFVAL_array(std::vector<unsigned char>& out, const Huff_Table& huff_table)
    {
        assert("huffman table should not have symbols above 255" && huff_table.size()<=256);

        // Can assume Huffman table is cannonical
        // Means for a given code size symbols appear with increasing codes,
        // for any i<j, table[i]<table[j]
        // So just need to apply the symbols in the order they appear, one code
        // size at a time
        const unsigned char min_code_size{ 1 };
        const unsigned char max_code_size{ 16 };

        for (size_t code_size=min_code_size; code_size<=max_code_size; code_size++)
        {
            for (size_t bs_ind = 0; bs_ind < huff_table.size(); bs_ind++)
            {
                const auto& bs{ huff_table[bs_ind] };

                if (bs.size()==code_size)
                {
                    // The symbol corresponds to index of the Huffman table
                    out.push_back(static_cast<unsigned char>(bs_ind));
                }
            }
        }
    }

    void append_huff_table_data(std::vector<unsigned char>& out, const Huff_Table_Ref& huff_table)
    {
        // first append the table class and destination, table class is most significant
        // table_class should be 0 for DC, 1 for AC
        unsigned int table_class{ huff_table.type==Huff_Table_Ref::Huff_Table_Type::AC };

        append_composite_byte(out, table_class, huff_table.destination_ind);

        // Append the BITS list
        append_BITS_array(out, huff_table.table);

        // Finaly append the HUFFVAL list
        append_HUFFVAL_array(out, huff_table.table);
    }

    void append_huff_table_marker_segment(std::vector<unsigned char>& out, std::vector<Huff_Table_Ref> tables)
    {
        // Append the DHT marker
        append_marker(out, Marker::Define_Huffman_Table);

        // Next comes the length of the table segment but we don't know that yet
        // So record the index of where the length should go and push two bytes
        const size_t length_ind{ out.size() };
        out.resize(out.size() + 2);

        // Now add each table in turn
        for (auto& table : tables)
        {
            append_huff_table_data(out, table);
        }
        
        // Finally go back and fill in the length, most significant byte first
        size_t bytes_appended{ out.size()-length_ind };

        out[length_ind] = static_cast<unsigned char>(bytes_appended >> 8);
        out[length_ind+1] = static_cast<unsigned char>(bytes_appended & 0xFF);
    }

    void append_scan_header(std::vector<unsigned char>& out, const std::vector<Comp_Info>& comp_infos)
    {
        // Append Start Of Scan marker
        append_marker(out, Marker::Start_Of_Scan);

        // Determine length
        unsigned int scan_header_length{ 6 + 2 * static_cast<unsigned int>(comp_infos.size()) };
        append_two_bytes(out, scan_header_length);

        // Add number of components in the scan
        out.push_back(static_cast<unsigned char>(comp_infos.size()));

        // Add component specific data
        for (size_t comp_ind = 0; comp_ind < comp_infos.size(); comp_ind++)
        {
            // Index of component in scan
            out.push_back(static_cast<unsigned char>(comp_ind));

            // Composite byte of DC and AC table destinations
            // DC index is most significant
            size_t dc_destination{ comp_infos[comp_ind].DC_Huff_table_ind };
            size_t ac_destination{ comp_infos[comp_ind].AC_Huff_table_ind };

            append_composite_byte(out, dc_destination, ac_destination);
        }

        // Start of spectral or predictor selection
        // Set to zero for sequential DCT
        out.push_back(0);

        // End of spectral selection
        // Set to 63 for sequential DCT
        out.push_back(63);
        
        // Composite byte of successive approximation bit position high & low
        // Both high and low are set to zero in sequential DCT
        out.push_back(0x00);
    }

    void append_frame_header(std::vector<unsigned char>& out, unsigned int Y, unsigned int X, 
        const std::vector<Comp_Info>& comp_infos)
    {
        // Append Start Of Frame 0 marker, corresponds to baseline DCT
        append_marker(out, Marker::Start_Of_Frame_0_Baseline_DCT);

        // Length of frame header
        unsigned int frame_header_length{ 8 + 3 * static_cast<unsigned int>(comp_infos.size()) };
        
        append_two_bytes(out, frame_header_length);

        // Sample precision in bits, 1 byte
        out.push_back(8);

        // Height and width of image, 2 bytes each
        append_two_bytes(out, Y);

        append_two_bytes(out, X);

        // Number of components in frame, 1 byte
        out.push_back(static_cast<unsigned char>(comp_infos.size()));

        // Component specific stuff
        for (size_t comp_ind = 0; comp_ind < comp_infos.size(); comp_ind++)
        {
            // Component identifier, 1 byte
            out.push_back(static_cast<unsigned char>(comp_ind));

            // Composite byte containing H sampling factor (most significant) and V sampling
            // factor (least significant)
            unsigned int H{ comp_infos[comp_ind].H }, V{ comp_infos[comp_ind].V };
            append_composite_byte(out, H, V);

            // Quantization table destination index, 1 byte
            out.push_back(static_cast<unsigned char>(comp_infos[comp_ind].q_table_ind)); 
        }
    }

    void encode_scan(std::vector<unsigned char>& out, const std::vector<Coefficient>& encoded_coeffs, 
        const std::vector<Comp_Info>& comp_infos, const std::vector<Huff_Table>& dc_tables, const std::vector<Huff_Table>& ac_tables)
    {
        // Append scan header
        append_scan_header(out, comp_infos);

        // Now need encode the bit string
        Bit_String bs;

        for (auto &coeff : encoded_coeffs)
        {
            size_t comp_ind{ coeff.comp_ind };

            if (coeff.type==Coefficient_Type::DC)
            {
                const Huff_Table& table{ dc_tables[comp_infos[comp_ind].DC_Huff_table_ind] };
                bs.extend(table[coeff.RS]);
            }
            else
            {
                const Huff_Table& table{ ac_tables[comp_infos[comp_ind].AC_Huff_table_ind] };
                bs.extend(table[coeff.RS]);
            }
            
            unsigned int ssss{ 0xFU & coeff.RS };
            
            // Note the spec says diff should be represented in 12 bit two's-complement.
            // The two's-complement of a negative number can therefore be computed as 
            // (~n) + 1. Note the two's-complement of 0 is 0.
            // So for the cas of diff<0 we can just append last ssss bits of ~diff.
            // Since diff=0 corresponds to ssss=0, we don't need to worry about it
            if (coeff.value>0)
            {
                bs.append_last_ssss_bits(coeff.value, ssss);
            }
            else
            {
                bs.append_last_ssss_bits(~static_cast<unsigned int>(-coeff.value), ssss);
            }
        }

        // Now need to add padding if necessary and then perform byte stuffing
        while (bs.size() % CHAR_BIT != 0)
        {
            bs.append_bit(1);
        }

        for (auto c_ptr=bs.begin_bytes(); c_ptr!=bs.end_bytes(); ++c_ptr)
        {
            out.push_back(*c_ptr);

            // If a 0xFF byte occurs, we must append a 0x00 byte after it
            if (*c_ptr==0xFF)
            {
                out.push_back(0x00);
            }
        }
    }

    void add_coeffs_stats(
        std::vector<Coefficient_Stats>& stats, const std::vector<Coefficient>& coeffs, const std::vector<Comp_Info>& comp_infos
    )
    {
        for (auto &coeff : coeffs)
        {
            // 
            size_t stat_ind;

            if (coeff.type==Coefficient_Type::DC)
            {
                stat_ind = comp_infos[coeff.comp_ind].DC_Huff_table_ind;
                stats[stat_ind].dc_stats[coeff.RS]++;
            }
            else
            {
                stat_ind = comp_infos[coeff.comp_ind].AC_Huff_table_ind;
                stats[stat_ind].ac_stats[coeff.RS]++;
            }
        }
    }

    void encode_frame(std::vector<unsigned char>& out, unsigned int Y, unsigned int X, std::vector<DU_Array<double>>& arrays, 
        const std::vector<Comp_Info>& comp_infos, const std::vector<Huff_Table>& provided_dc_tables, const std::vector<Huff_Table>& provided_ac_tables,
        const std::vector<Q_Table>& q_tables, bool optimize_huff)
    {
        // Encode the image data into the Coefficient representation
        // This is so we can optionally produce optimized Huffman tables for this particular image
        std::vector<Coefficient> encoded_coeffs{ encode_coeff_rep_sequential(arrays, comp_infos) };

        std::vector<Huff_Table> optimized_ac_tables, optimized_dc_tables;

        // Generate optimized Huffman tables if thet's what we want to do
        if (optimize_huff)
        {
            // First determine How many DC and AC tables we need to compute
            size_t max_dc_table_ind{ 
                std::max_element(comp_infos.begin(), comp_infos.end(), 
                    [](const Comp_Info& a, const Comp_Info& b)
                    {
                        return a.DC_Huff_table_ind<b.DC_Huff_table_ind;
                    }
                )->DC_Huff_table_ind
            };
            size_t max_ac_table_ind{ 
                std::max_element(comp_infos.begin(), comp_infos.end(), 
                    [](const Comp_Info& a, const Comp_Info& b)
                    {
                        return a.AC_Huff_table_ind<b.AC_Huff_table_ind;
                    }
                )->AC_Huff_table_ind
            };

            // Compute the RS statistics
            std::vector<Coefficient_Stats> stats(1+std::max(max_dc_table_ind, max_ac_table_ind));

            add_coeffs_stats(stats, encoded_coeffs, comp_infos);

            // Finnaly add the tables
            for (auto &coeff_stat : stats)
            {
                // generate DC table
                Huff_Table&& dc_table{ Huff_Table::gen_table_from_stats(coeff_stat.dc_stats) };
                optimized_dc_tables.push_back(std::move(dc_table));
                
                // generate AC table
                Huff_Table&& ac_table{ Huff_Table::gen_table_from_stats(coeff_stat.ac_stats) };
                optimized_ac_tables.push_back(std::move(ac_table));
            }
        }

        // Now set the table refs that point to which set of tables we should use
        auto& dc_tables{ optimize_huff ? optimized_dc_tables : provided_dc_tables };
        auto& ac_tables{ optimize_huff ? optimized_ac_tables : provided_ac_tables };

        // Append various tables starting with quantization tables       
        append_q_table_marker_segment(out, q_tables);

        // Now add Huffman tables
        std::vector<Huff_Table_Ref> huff_refs;

        for (auto &comp_info : comp_infos)
        {
            huff_refs.emplace_back(
                dc_tables[comp_info.DC_Huff_table_ind], 
                Huff_Table_Ref::Huff_Table_Type::DC,
                comp_info.DC_Huff_table_ind
            );
            huff_refs.emplace_back(
                ac_tables[comp_info.AC_Huff_table_ind], 
                Huff_Table_Ref::Huff_Table_Type::AC,
                comp_info.AC_Huff_table_ind
            );
        }
        
        append_huff_table_marker_segment(out, huff_refs);

        // Now append the frame header
        append_frame_header(out, Y, X, comp_infos);

        // Finally encode the scan
        encode_scan(out, encoded_coeffs, comp_infos, dc_tables, ac_tables);
    }

    DU_Array<double> convert_to_DU_Array(const Array_2d<double>& array_2d, const Comp_Info& comp_info)
    {   
        // Number of data units
        const size_t N_du{ array_2d.size() / du_size };
        // Width of array_2d in units of data untis
        const size_t width_du{ array_2d.shape()[1] / du_width };

        const size_t H{ comp_info.H }, V{ comp_info.V };

        DU_Array<double> du_array(N_du);

        for (size_t ind = 0; ind < N_du; ind++)
        {
            // Index of the MCU du_array[ind] exists in
            const size_t mcu_ind{ ind / (H*V) };
            const size_t offset{ ind % (H*V) };

            // position of first data unit in the mcu in array_2d in units of data units 
            const size_t mcu_i{ ((mcu_ind*H) / width_du) * V };
            const size_t mcu_j{ (mcu_ind*H) % width_du };

            // position of the data unit we're actually interested in,
            // in units of data units
            const size_t du_i{ mcu_i + offset / H };
            const size_t du_j{ mcu_j + offset % H };

            // Indices in array_2d of the data unit we're copying
            const size_t i_begin{ du_height*du_i };
            const size_t j_begin{ du_width*du_j };

            for (size_t i = 0; i < du_height; i++)
            {
                for (size_t j = 0; j < du_width; j++)
                {
                    du_array(ind, i, j) = array_2d(i_begin+i, j_begin+j);
                }
            }
        }

        return du_array;
    }

    std::vector<unsigned char> encode_image(unsigned int Y, unsigned int X, const std::vector<Array_2d<double>>& arrays, 
        const std::vector<Comp_Info>& comp_infos, const std::vector<Huff_Table>& dc_tables, const std::vector<Huff_Table>& ac_tables,
        const std::vector<Q_Table>& q_tables, bool optimize_huff)
    {
        // Reshape the arrays into DU_Arrays
        std::vector<DU_Array<double>> du_arrays;

        for (size_t comp_ind = 0; comp_ind < comp_infos.size(); comp_ind++)
        {
            auto& array{ arrays[comp_ind] };
            auto& comp_info{ comp_infos[comp_ind] };

            du_arrays.push_back(std::move(convert_to_DU_Array(array, comp_info)));
        }
        
        // Perform numerical operations
        for (size_t comp_ind = 0; comp_ind < du_arrays.size(); comp_ind++)
        {
            auto& array{ du_arrays[comp_ind] };

            apply_level_shift(array);
            apply_DCT(array);

            auto& q_table{ q_tables[comp_infos[comp_ind].q_table_ind] };

            apply_quantization(array, q_table);
        }

        // Now we can encode the image
        // First apply the start of image marker, FFD8
        std::vector<unsigned char> encoded_image;

        append_marker(encoded_image, Marker::Start_Of_Image);

        // Append the encoded frame
        encode_frame(encoded_image, Y, X, du_arrays, comp_infos, dc_tables, ac_tables, q_tables, optimize_huff);

        // Lastly append the end of image marker, FFD9
        append_marker(encoded_image, Marker::End_Of_Image);

        return encoded_image;
    }

    Array_2d<double> enlarge_component(const Array_2d<double>& orig_comp, unsigned int V, unsigned int H)
    {
        const size_t cur_height{ orig_comp.shape()[0] }, cur_width{ orig_comp.shape()[1] };
        // Size of MCU in units of samples
        const size_t mcu_height{ du_height*V }, mcu_width{ du_width*H };

        // The height/width of the enlarged component
        const size_t required_height{ cur_height % mcu_height ? (1+cur_height/mcu_height)*mcu_height : cur_height };
        const size_t required_width{ cur_width % mcu_width ? (1+cur_width/mcu_width)*mcu_width : cur_width };
        
        // Both height and width requirements met, just return a copy
        if (cur_height==required_height && cur_width==required_width)
        {
            return orig_comp;
        }

        // Need to enlarge
        Array_2d<double> enlarged_component{ required_height, required_width };

        for (size_t i = 0; i < required_height; i++)
        {
            for (size_t j = 0; j < required_width; j++)
            {
                if (i<cur_height && j<cur_width)
                {
                    // Copy the original data
                    enlarged_component(i, j) = orig_comp(i, j);
                }
                else if (i>=cur_height && j<cur_width)
                {
                    // Rows beneath the original data
                    enlarged_component(i, j) = orig_comp(cur_height-1, j);
                }
                else if (i<cur_height && j>=cur_width)
                {
                    // Columns to the right of the original data
                    enlarged_component(i, j) = orig_comp(i, cur_width-1);
                }
                else
                {
                    // Rectangle in the far right, set to 128 as that seems a neutral choice
                    enlarged_component(i, j) = 128;
                }
            }
        }
        
        return enlarged_component;
    }

    std::tuple<Array_2d<double>&, Array_2d<double>&, Array_2d<double>&> colour_transform(
        Array_2d<double>& red, Array_2d<double>& green, Array_2d<double>& blue
    )
    {
        assert("red and green components were not the same size" && red.size()==green.size());
        assert("green and blue components were not the same size" && green.size()==blue.size());

        auto& Y_arr{ red };
        auto& U_arr{ green };
        auto& V_arr{ blue };

        // Note round() is defined as floor(x+0.5) in this context
        auto round = [](double n){
            return std::floor(n+0.5);
        };

        constexpr size_t N_colours{ 3 };

        std::array<std::array<double, N_colours>, N_colours> colour_mat{{
            {0.299, 0.587, 0.114},
            {-0.299 / 1.772, -0.587 / 1.772,  0.886 / 1.772},
            { 0.701 / 1.402, -0.587 / 1.402, -0.114 / 1.402}
        }};

        std::array<double, N_colours> biases{ 0, 128, 128 };

        std::array<double, N_colours> RGB_cols, YUV_cols;

        for (size_t ind = 0; ind < red.size(); ind++)
        {
            RGB_cols = {red[ind], green[ind], blue[ind]};
            
            for (size_t yuv_ind = 0; yuv_ind < N_colours; yuv_ind++)
            {
                double& cur_YUV{ YUV_cols[yuv_ind] };

                cur_YUV = biases[yuv_ind];

                for (size_t rgb_ind = 0; rgb_ind < N_colours; rgb_ind++)
                {
                    cur_YUV += colour_mat[yuv_ind][rgb_ind]*RGB_cols[rgb_ind];
                }

                cur_YUV = std::min(std::max(0.0, round(cur_YUV)), 255.0);
            }
            
            Y_arr[ind] = YUV_cols[0];
            U_arr[ind] = YUV_cols[1];
            V_arr[ind] = YUV_cols[2];
        }

        return {Y_arr, U_arr, V_arr};
    }

    void subsample_component_4_2_0(Array_2d<double>& component)
    {
        size_t height{ component.shape()[0] }, width{ component.shape()[1] };

        // In 4:2:0 subsampling we average over 2x2 data unit sized blocks
        // Check height and width are both multiples of 16
        if (height % (2*du_height))
        {
            throw std::invalid_argument("Component height must be a multiple of 16 for 4:2:0 subsampling");
        }
        else if (width % (2*du_width))
        {
            throw std::invalid_argument("Component width must be a multiple of 16 for 4:2:0 subsampling");
        }

        // Perform the subsampling
        Array_2d<double> component_ss{ height / 2, width / 2 };

        for (size_t i = 0; i < component_ss.shape()[0]; i++)
        {
            for (size_t j = 0; j < component_ss.shape()[1]; j++)
            {
                component_ss(i, j) = (
                    component(2*i,   2*j  ) +
                    component(2*i,   2*j+1) +
                    component(2*i+1, 2*j  ) +
                    component(2*i+1, 2*j+1)
                ) / 4.0;
            }
        }

        component = std::move(component_ss);
    }

    void subsample_component_4_2_2(Array_2d<double>& component)
    {
        size_t height{ component.shape()[0] }, width{ component.shape()[1] };

        // In 4:2:0 subsampling we average over 1x2 data unit sized blocks, that is
        // 1 pixel tall, 2 pixels wide
        // Check height is multiple of 8, width mutliple of 16
        if (height % du_height)
        {
            throw std::invalid_argument("Component height must be a multiple of 8 for 4:2:2 subsampling");
        }
        else if (width % (2*du_width))
        {
            throw std::invalid_argument("Component width must be a multiple of 16 for 4:2:2 subsampling");
        }

        // Perform the subsampling
        Array_2d<double> component_ss{ height, width / 2 };

        for (size_t i = 0; i < component_ss.shape()[0]; i++)
        {
            for (size_t j = 0; j < component_ss.shape()[1]; j++)
            {
                component_ss(i, j) = (
                    component(i,   2*j  ) +
                    component(i,   2*j+1)
                ) / 2.0;
            }
        }

        component = std::move(component_ss);
    }

    void subsample_component(Array_2d<double>& component, Subsampling ss)
    {
        switch (ss)
        {
        case Subsampling::ss_4_4_4:
            // 4:4:4 subsampling leaves component untouched
            break;
        case Subsampling::ss_4_2_2:
            subsample_component_4_2_2(component);
            break;
        case Subsampling::ss_4_2_0:
            subsample_component_4_2_0(component);
            break;
        default:
            throw std::invalid_argument("Unsupported subsampling");
        }
    }
}

