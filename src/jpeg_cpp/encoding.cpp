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
        const std::array<size_t, 64> zz{ 
            0, 1, 8, 16, 9, 2, 3, 10, 17, 24, 32, 
            25, 18, 11, 4, 5, 12, 19, 26, 33, 40, 
            48, 41, 34, 27, 20, 13, 6, 7, 14, 21, 
            28, 35, 42, 49, 56, 57, 50, 43, 36, 
            29, 22, 15, 23, 30, 37, 44, 51, 58, 
            59, 52, 45, 38, 31, 39, 46, 53, 60, 
            61, 54, 47, 55, 62, 63
         };
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
}
