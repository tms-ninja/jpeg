#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "jpeg_cpp/array.h"
#include "jpeg_cpp/encoding.h"

using namespace Catch::Matchers;
using namespace JPEG;

// Generates a DU_Array containing 2 data units with numbers from 0 to 127
DU_Array<double> gen_arange_DU_Array()
{
    DU_Array<double> array{ 2 };

    for (size_t ind = 0; ind < array.size(); ind++)
    {
        array[ind] = ind;
    }
    
    return array;
}

TEST_CASE( "apply_level_shift()::level shift values in DU_Array", "[apply_level_shift()]" ) {
    DU_Array<double> array{ gen_arange_DU_Array() };
    
    apply_level_shift(array);

    // Verify all elements have been reduced by 128
    for (size_t ind = 0; ind < array.size(); ind++)
    {
        double expected{ static_cast<double>(ind) - 128.0 };
        
        REQUIRE_THAT( array[ind], WithinRel(expected) );
    }
}


