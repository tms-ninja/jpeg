#include "jpeg_cpp/coefficient.h"

namespace JPEG
{
    bool operator==(const Coefficient& coeff1, const Coefficient& coeff2)
    {
        return (
            coeff1.value == coeff2.value &&
            coeff1.RS == coeff2.RS &&
            coeff1.type == coeff2.type &&
            coeff1.comp_ind == coeff2.comp_ind
        );
    }

    bool operator!=(const Coefficient& coeff1, const Coefficient& coeff2)
    {
        return !(coeff1==coeff2);
    }

    std::ostream& operator<<(std::ostream& out, const Coefficient& coeff)
    {
        out << "{value=" << coeff.value << ", RS=" << coeff.RS << ", type=";

        if (coeff.type==Coefficient_Type::DC)
        {
            out << "DC";
        }
        else
        {
            out << "AC";
        }

        out << ", comp_ind=" << coeff.comp_ind << '}';

        return out;
    }
}
