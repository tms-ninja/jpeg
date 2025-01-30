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
        std::array<double, 64> double_buffer;
        const auto arr_shape{ array.shape() };

        for (size_t du_ind = 0; du_ind < arr_shape[0]; du_ind++)
        {
            // Compute data unit*DCT matrix
            mat_mul(&array[du_ind*64], DCT_matrix.data(), &double_buffer[0], arr_shape[1]);

            // Compute transpose(DCT matrix) * termporary
            mat_mul(DCT_matrix_transpose.data(), &double_buffer[0], &array[du_ind*64], arr_shape[1]);
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
                q_table_corrected(i, j) = q_table(i, j);

                if (i==0)
                {
                    q_table_corrected(i, j) *= std::sqrt(2.0);
                }
                
                if (j==0)
                {
                    q_table_corrected(i, j) *= std::sqrt(2.0);
                }
            }
        }

        // Now we can apply quantization
        for (size_t du_ind = 0; du_ind < array.shape()[0]; du_ind++)
        {
            for (size_t i = 0; i < array.shape()[1]; i++)
            {
                for (size_t j = 0; j < array.shape()[2]; j++)
                {
                    array(du_ind, i, j) = std::round(array(du_ind, i, j) / q_table_corrected(i, j));
                }
            }
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

    void encode_DC_coeff(Bit_String& bs, int diff, const Huff_Table& huff_table)
    {
        unsigned int ssss{ compute_ssss(std::abs(diff)) };
        
        // Append the Huffman code corresponding to the ssss value
        bs.extend(huff_table[ssss]);

        // Now need to append last ssss bits of:
        // - diff if diff>0
        // - (diff-1) if diff<0
        // Note the spec says diff should be represented in 12 bit two's-complement.
        // The two's-complement of a negative number can therefore be computed as 
        // (~n) + 1. Note the two's-complement of 0 is 0.
        // So for the cas of diff<0 we can just append last ssss bits of ~diff.
        // Since diff=0 corresponds to ssss=0, we don't need to worry about it
        if (diff==0)
        {
            return;
        }

        if (diff>0)
        {
            bs.append_last_ssss_bits(diff, ssss);
        }
        else
        {
            bs.append_last_ssss_bits(~static_cast<unsigned int>(-diff), ssss);
        }
    }
    
    void encode_AC_coeff(Bit_String& bs, int coeff, unsigned int rrrr, const Huff_Table& huff_table)
    {
        assert("coeff should not be zero" && coeff!=0);
        assert("rrrr must be less than 16" && rrrr<16);

        unsigned int ssss{ compute_ssss(std::abs(coeff)) };
        unsigned int rrrrssss{ (rrrr << 4) + ssss };

        // Append the Huffman code encoding rrrrssss to the output
        bs.extend(huff_table[rrrrssss]);

        // Like with the DC coefficient we need to simply append last ssss bits of ac_coeff
        // if ac_coeff>0 or ~(-ac_coeff) if it's negative
        unsigned int value_to_encode{ coeff>0 ? coeff : ~static_cast<unsigned int>(-coeff) };

        bs.append_last_ssss_bits(value_to_encode, ssss);
    }

    void encode_AC_coeffs(Bit_String& bs, const DU_Array<double>& du_array, size_t du_ind,
                            const Huff_Table& huff_table)
    {
        // Offsets for the zig-zag pattern
        const auto& zz{ zig_zag_order };

        // Pointer to the start of the data unit we are encoding
        const double* du_ptr{ du_array.data() + du_ind*64 };
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
                    bs.extend(huff_table[end_of_block]);
                }

                r++;
                continue;
            }

            // Note r could be much larger than 15 and might need to encode more 
            // than one ZRL
            while (r>15)
            {
                // Encode ZRL
                bs.extend(huff_table[zrl]);
                r -= 16;
            }

            // Now encode the current, non-zero coefficient
            encode_AC_coeff(bs, value_to_encode, static_cast<unsigned char>(r), huff_table);

            // Reset the run length of zeros
            r = 0;
        }
    }

    int encode_data_unit_sequential(Bit_String& bs, const DU_Array<double>& du_array, size_t du_ind, int prev_dc,
                        const Huff_Table& huff_table_dc, const Huff_Table& huff_table_ac)
    {
        // First encode the DC coefficient
        int dc_coeff{ static_cast<int>(du_array(du_ind, 0, 0)) };
        // Its the difference between the DC coeff in this unit and the previous one that
        // is encoded
        int dc_diff{ dc_coeff-prev_dc };

        encode_DC_coeff(bs, dc_diff, huff_table_dc);

        // Now encode the AC coefficients
        encode_AC_coeffs(bs, du_array, du_ind, huff_table_ac);

        return dc_coeff;
    }

    void append_q_table_marker_segment(std::vector<unsigned char>& out, const std::vector<Q_Table>& q_tables, std::vector<unsigned int>& destination_indices)
    {
        assert("Number of quantization tables and destination indices must be the same" && q_tables.size()==destination_indices.size());

        // Append quantization table marker, 0xFFDB, most significant byte first
        out.push_back(0xFF);
        out.push_back(0xDB);

        // Now determine the length of the marker segment
        size_t length{ 2 + 65 * q_tables.size() };

        append_two_bytes(out, length);
        // out.push_back(length >> 8);
        // out.push_back(0xFF & length);

        // Now append each table in turn
        for (size_t table_ind = 0; table_ind < q_tables.size(); table_ind++)
        {
            const Q_Table& table{ q_tables[table_ind] };

            // Table precision, 0 for 8 bit images
            unsigned int Pq{ 0 };

            // Table destination identifier
            unsigned int Tq{ destination_indices[table_ind] };
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

    void append_mcu(Bit_String& bs, std::vector<int>& prev_dc, std::vector<size_t>& du_ind, const std::vector<DU_Array<double>>& arrays, 
        const std::vector<Comp_Info>& comp_infos, const std::vector<Huff_Table>& dc_tables, const std::vector<Huff_Table>& ac_tables)
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
                    bs, 
                    data_units, 
                    du_ind[comp_ind],
                    prev_dc[comp_ind],
                    dc_tables[comp_infos[comp_ind].DC_Huff_table_ind],
                    ac_tables[comp_infos[comp_ind].AC_Huff_table_ind]
                );

                // Move on to next data unit
                du_ind[comp_ind]++;
            }
        }
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
        assert("huffman table should not have symbols above 255" && huff_table.size()<255);

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
        out.push_back(0xFF);
        out.push_back(0xC4);

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
        out.push_back(0xFF);
        out.push_back(0xDA);

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
        out.push_back(0xFF);
        out.push_back(0xC0);

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

    void encode_scan(std::vector<unsigned char>& out, std::vector<DU_Array<double>>& arrays, 
        const std::vector<Comp_Info>& comp_infos, const std::vector<Huff_Table>& dc_tables, const std::vector<Huff_Table>& ac_tables)
    {
        // Append scan header
        append_scan_header(out, comp_infos);
        
        // Entropy encode each MCU in turn
        // Note previous DC should be initialized to 0
        Bit_String bs;
        std::vector<int> prev_dc(arrays.size());
        std::vector<size_t> du_ind(arrays.size());

        while (du_ind[0]<arrays[0].shape()[0])
        {
            append_mcu(bs, prev_dc, du_ind, arrays, comp_infos, dc_tables, ac_tables);
        }

        // Now need to add padding if necessary and then perform byte stuffing
        while (bs.size() % 8 != 0)
        {
            bs.append_bit(1);
        }

        const auto end{ bs.data() + bs.size() / bits_per_char };

        for (auto c_ptr=bs.data(); c_ptr!=end; ++c_ptr)
        {
            out.push_back(*c_ptr);

            // If a 0xFF byte occurs, we must append a 0x00 byte after it
            if (*c_ptr==0xFF)
            {
                out.push_back(0x00);
            }
        }
    }

    void encode_frame(std::vector<unsigned char>& out, unsigned int Y, unsigned int X, std::vector<DU_Array<double>>& arrays, 
        const std::vector<Comp_Info>& comp_infos, const std::vector<Huff_Table>& dc_tables, const std::vector<Huff_Table>& ac_tables,
        const std::vector<Q_Table>& q_tables)
    {
        // Append various tables starting with quantization tables
        std::vector<unsigned int> q_table_inds;

        for (auto &comp_info : comp_infos)
        {
            q_table_inds.push_back(static_cast<unsigned char>(comp_info.q_table_ind));
        }
        
        append_q_table_marker_segment(out, q_tables, q_table_inds);

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
        encode_scan(out, arrays, comp_infos, dc_tables, ac_tables);
    }

    DU_Array<double> convert_to_DU_Array(const Array_2d<double>& array_2d, const Comp_Info& comp_info)
    {   
        // Number of data units
        const size_t N_du{ array_2d.size() / 64 };
        // Width of array_2d in units of data untis
        const size_t width_du{ array_2d.shape()[1] / 8 };

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
            const size_t i_begin{ 8*du_i };
            const size_t j_begin{ 8*du_j };

            for (size_t i = 0; i < 8; i++)
            {
                for (size_t j = 0; j < 8; j++)
                {
                    du_array(ind, i, j) = array_2d(i_begin+i, j_begin+j);
                }
            }
        }

        return du_array;
    }

    std::vector<unsigned char> encode_image(unsigned int Y, unsigned int X, const std::vector<Array_2d<double>>& arrays, 
        const std::vector<Comp_Info>& comp_infos, const std::vector<Huff_Table>& dc_tables, const std::vector<Huff_Table>& ac_tables,
        const std::vector<Q_Table>& q_tables)
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

        encoded_image.push_back(0xFF);
        encoded_image.push_back(0xD8);

        // Append the encoded frame
        encode_frame(encoded_image, Y, X, du_arrays, comp_infos, dc_tables, ac_tables, q_tables);

        // Lastly append the end of image marker, FFD9
        encoded_image.push_back(0xFF);
        encoded_image.push_back(0xD9);

        return encoded_image;
    }

    Array_2d<double> enlarge_component(const Array_2d<double> orig_comp, unsigned int V, unsigned int H)
    {
        const size_t cur_height{ orig_comp.shape()[0] }, cur_width{ orig_comp.shape()[1] };
        const size_t mcu_height{ 8*V }, mcu_width{ 8*H };

        // The height/width of the enlarged component
        const size_t required_height{ cur_height % mcu_height ? (1+cur_height/mcu_height)*mcu_height : cur_height };
        const size_t required_width{ cur_width % mcu_width ? (1+cur_width/mcu_width)*mcu_width : cur_width };
        
        // Both height and width requirements met, just return a copy
        if (cur_height==required_height && cur_width==required_width)
        {
            return Array_2d<double>{orig_comp};
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
}

