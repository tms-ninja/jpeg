#ifndef BIT_STRING_H
#define BIT_STRING_H

#include <climits>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace JPEG {

    // Check unsigned char have exactly 8 bits
    static_assert(CHAR_BIT==8, "char type has more than 8 bits!");

    // bit selector for the most significant bit in an 8 bit byte
    constexpr unsigned int MSB_selector{ 0b10'000'000u };


    /// @brief Proxy class to allow subscripting bits in Bit_String
    class Bit_Ref
    {
    private:
        unsigned char& char_ref;
        const size_t bit_offset;
        const size_t bit_selector;

    public:
        Bit_Ref(unsigned char& c, size_t offset);
        
        Bit_Ref& operator=(bool bit);
        operator bool();
    };


    /// @brief Represents a String of Bits
    /// Note first bit in an 8 bit byte is the most significant bit
    /// This should be the correct order for the JPEG standard
    class Bit_String {
    private:
        std::vector<unsigned char> buffer;
        size_t buffer_size;

    public:
        /// @brief Creates an empty Bit_String
        Bit_String();

        /// @brief Creates a zeroed Bit_String of the given size
        /// @param size Number of bits of the new Bit_String
        Bit_String(size_t size);

        /// @brief Iterator to the first byte in the Bit_String. If the Bit_String does not represent an
        /// integer number of bytes, the least significant bits of the last byte are padded with zeros.
        /// @return Iterator to the first byte in the Bit_String
        unsigned char* begin_bytes() { return buffer.data(); }

        /// @brief Iterator to one past the last byte in the Bit_String
        /// @return Iterator to one past the last byte in the Bit_String
        unsigned char* end_bytes() { return buffer.data()+buffer.size(); }

        /// @brief Creates a Bit_String from a string
        /// @param s String containing '0' and '1'
        /// @return New Bit_String
        Bit_String(const std::string& s);

        /// @brief Number of bits in the Bit_String
        /// @return Number of bits
        size_t size() const;

        /// @brief Pointer to the internal buffer containing the bits
        /// @return Pointer to the internal buffer containing the bits
        const unsigned char* data() const;

        // Subscript operators
        Bit_Ref operator[](size_t ind);
        bool operator[](size_t ind) const;

        // Comparison operators
        friend bool operator==(const Bit_String& bs1, const Bit_String& bs2);
        friend bool operator!=(const Bit_String& bs1, const Bit_String& bs2);

        // I/O output operator
        friend std::ostream& operator<< (std::ostream& out, const Bit_String& bs);

        /// @brief Appends a bit
        /// @param bit Bit to append 
        void append_bit(bool bit);

        /// @brief Appends last ssss bits of n, most significant bit first
        /// @param n Source of bits to append
        /// @param ssss Number of bits to append
        void append_last_ssss_bits(unsigned int n, unsigned int ssss);

        /// @brief Appends the given byte, most significant bit first
        /// @param byte Byte to append
        void append_byte(unsigned char byte);

        /// @brief Extends the given Bit_string with the contents of another Bit_String
        /// @param other_bs The other Bit_String from which bits are sourced
        void extend(const Bit_String& other_bs);
    };
}

#endif
