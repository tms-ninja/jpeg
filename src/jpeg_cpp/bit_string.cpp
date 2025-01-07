#include "jpeg_cpp/bit_string.h"

Bit_String::Bit_String() 
 : buffer_size{ 0 }
{

}

Bit_String::Bit_String(size_t size)
 : buffer_size{ size }
{
    buffer.resize(size/bits_per_char + (size % bits_per_char > 0));
}

size_t Bit_String::size() const
{
    return buffer_size;
}

const unsigned char* Bit_String::data() const
{
    return buffer.data();
}

Bit_Ref Bit_String::operator[](size_t ind)
{
    size_t byte_offset{ ind / bits_per_char }, bit_offset{ ind % bits_per_char };

    return Bit_Ref(buffer[byte_offset], bit_offset);
}

bool Bit_String::operator[](size_t ind) const
{
    size_t byte_offset{ ind / bits_per_char }, bit_offset{ ind % bits_per_char };
    size_t bit_selector{ MSB_selector >> bit_offset };

    return (buffer[byte_offset] & bit_selector);// >> (bits_per_char-bit_offset-1);
}

bool operator==(const Bit_String& bs1, const Bit_String& bs2)
{
    // Check if they're actually the same object
    if (&bs1==&bs2)
    {
        return true;
    }

    // Check if they have different sizes
    if (bs1.size() != bs2.size())
    {
        return false;
    }

    // Now check each elements
    for (size_t ind=0; ind<bs1.size(); ++ind)
    {
        if (bs1[ind] != bs2[ind])
        {
            return false;
        }
    }

    return true;
}

bool operator!=(const Bit_String& bs1, const Bit_String& bs2)
{
    return !(bs1==bs2);
}

std::ostream& operator<< (std::ostream& out, const Bit_String& bs)
{
    for (size_t ind=0; ind<bs.size(); ++ind)
    {
        out << static_cast<int>(bs[ind]);
    }

    return out;
}

void Bit_String::append_bit(bool bit)
{
    // check to see if we need to add a new unsigned char to the buffer
    if (buffer_size % bits_per_char == 0)
    {
        buffer.push_back(0u);
    }

    // Set the bit and update the size of the buffer
    (*this)[buffer_size] = bit;
    buffer_size++;
}

void Bit_String::append_last_ssss_bits(unsigned int n, unsigned int ssss)
{
    if (ssss==0)
    {
        return;
    }

    for (int bit_ind=ssss-1; bit_ind>=0; --bit_ind)
    {
        append_bit((n & (1u<<bit_ind)) >> bit_ind);
    }
}

void Bit_String::append_byte(unsigned char byte)
{
    // Check to see if we can just append the byte
    if (buffer_size % bits_per_char == 0)
    {
        buffer.push_back(byte);
    
        // Don't forget to increment the number of bits
        buffer_size += bits_per_char;
        return;
    }

    // Need to do one bit at a time
    for (size_t bit_ind = 0; bit_ind < bits_per_char; ++bit_ind)
    {
        append_bit(byte & (MSB_selector >> bit_ind));
    }
}

void Bit_String::extend(const Bit_String& other_bs)
{
    // If the size of this Bit_String is a multiple of 8, we can just extend the vector
    if (buffer_size % bits_per_char == 0)
    {
        for (auto byte: other_bs.buffer)
        {
            buffer.push_back(byte);
        }

        // Don't forget to update the size of this Bit_string
        buffer_size += other_bs.buffer_size;
        return;
    }

    // Deal with the case where buffer_size is not a multiple of 8
    // Need to do one bit at a time
    for (size_t bit_ind = 0; bit_ind < other_bs.size(); ++bit_ind)
    {
        append_bit(other_bs[bit_ind]);
    }
}

Bit_String Bit_String::from_string(const std::string& s)
{
    Bit_String bs{};

    for (auto c: s)
    {
        switch (c)
        {
        case '0':
            bs.append_bit(0u);
            break;    
        case '1':
            bs.append_bit(1u);
            break;  
        default:
            std::string what_msg{"'"};
            what_msg += c;
            what_msg += "' is not a valid character for a Bit_String (expected '0' or '1')";

            throw std::domain_error(what_msg);
        }
    }

    return bs;
}

