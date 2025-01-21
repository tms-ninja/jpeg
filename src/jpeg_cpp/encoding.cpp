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
            encode_AC_coeff(bs, value_to_encode, r, huff_table);

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

        out.push_back(length >> 8);
        out.push_back(0xFF & length);

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
            out.push_back((Pq << 4) + Tq);

            // Now need to append the value of the quantization table in zig-zag order
            const auto& zz{ zig_zag_order };

            for (size_t ind = 0; ind < table.size(); ind++)
            {   
                assert("Quantization table values must be less than 256" && table[zz[ind]]);
                out.push_back(table[zz[ind]]);
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

    void append_BITS_array(std::vector<unsigned char>& out, Huff_Table& huff_table)
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

    void append_HUFFVAL_array(std::vector<unsigned char>& out, Huff_Table& huff_table)
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
                    out.push_back(bs_ind);
                }
            }
        }
    }

    void append_huff_table_data(std::vector<unsigned char>& out, Huff_Table_Ref& huff_table)
    {
        // first append the table class and destination, table class is most significant
        // table_class should be 0 for DC, 1 for AC
        unsigned int table_class{ huff_table.type==Huff_Table_Ref::Huff_Table_Type::AC };

        out.push_back((table_class << 4) + huff_table.destination_ind);

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

        out[length_ind] = bytes_appended >> 8;
        out[length_ind+1] = bytes_appended & 0xFF;
    }

    void append_scan_header(std::vector<unsigned char>& out, const std::vector<Comp_Info>& comp_infos)
    {
        // Append Start Of Scan marker
        out.push_back(0xFF);
        out.push_back(0xDA);

        // Determine length
        unsigned int scan_header_length{ 6 + 2 * static_cast<unsigned int>(comp_infos.size()) };
        out.push_back(scan_header_length >> 8);
        out.push_back(scan_header_length & 0xFF);

        // Add number of components in the scan
        out.push_back(comp_infos.size());

        // Add component specific data
        for (size_t comp_ind = 0; comp_ind < comp_infos.size(); comp_ind++)
        {
            // Index of component in scan
            out.push_back(comp_ind);

            // Composite byte of DC and AC table destinations
            // DC index is most significant
            size_t dc_destination{ comp_infos[comp_ind].DC_Huff_table_ind };
            size_t ac_destination{ comp_infos[comp_ind].AC_Huff_table_ind };
            out.push_back((dc_destination << 4) + ac_destination);
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
        out.push_back(frame_header_length >> 8);
        out.push_back(frame_header_length & 0xFF);

        // Sample precision in bits, 1 byte
        out.push_back(8);

        // Height and width of image, 2 bytes each
        out.push_back(Y >> 8);
        out.push_back(Y & 0xFF);

        out.push_back(X >> 8);
        out.push_back(X & 0xFF);

        // Number of components in frame, 1 byte
        out.push_back(comp_infos.size());

        // Component specific stuff
        for (size_t comp_ind = 0; comp_ind < comp_infos.size(); comp_ind++)
        {
            // Component identifier, 1 byte
            out.push_back(comp_ind);

            // Composite byte containing H sampling factor (most significant) and V sampling
            // factor (least significant)
            unsigned int H{ comp_infos[comp_ind].H }, V{ comp_infos[comp_ind].V };
            out.push_back((H << 4) + V);

            // Quantization table destination index, 1 byte
            out.push_back(comp_infos[comp_ind].q_table_ind); 
        }
    }
}

