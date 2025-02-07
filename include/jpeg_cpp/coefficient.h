#ifndef COEFFICIENT_H
#define COEFFICIENT_H

#include <iostream>
#include <vector>

namespace JPEG
{
    enum class Coefficient_Type
    {
        DC,
        AC
    };
    
    /// @brief Represents a quantized coefficient
    struct Coefficient
    {
        // For 8 bit images the result of the quantization step is represented using a 
        // 12 bit integer, so short must be at least 2 bytes
        static_assert(sizeof(short)>=2, "Short must be at least 2 bytes");

        short value;
        unsigned char RS;
        Coefficient_Type type;
        unsigned char comp_ind;

        /// @brief Constructs a Coefficient
        /// @param value Value of the coefficient
        /// @param RS The RRRRSSSS byte representation of the coefficient
        /// @param type The coefficient type, DC or AC
        /// @param comp_ind The index of the component this coefficient corresponds to
        Coefficient(short value, unsigned char RS, Coefficient_Type type, unsigned char comp_ind)
        : value{ value }, RS{ RS }, type{ type }, comp_ind{ comp_ind }
        {

        }

        friend bool operator==(const Coefficient& coeff1, const Coefficient& coeff2);
        friend bool operator!=(const Coefficient& coeff1, const Coefficient& coeff2);

        friend std::ostream& operator<<(std::ostream& out, const Coefficient& coeff);
    };

    /// @brief Represents statistics about a component
    struct Coefficient_Stats
    {

        std::vector<unsigned int> dc_stats;
        std::vector<unsigned int> ac_stats;

        Coefficient_Stats() : dc_stats(256), ac_stats(256) {}
    };
}

#endif
