#include <catch2/catch_test_macros.hpp>
#include "jpeg_cpp/bit_string.h"
#include "jpeg_cpp/huff_table.h"

#include <string>
#include <vector>

using namespace JPEG;

TEST_CASE( "Huff_Table::construction of Huff_Tables", "[Huff_Table]" ) {
    size_t expected_size{ 5 };
    
    SECTION( "construction of given size" ) {
        Huff_Table huff_table_actual( expected_size );

        REQUIRE( huff_table_actual.size() == expected_size );

        // Check elements are default constructed
        for (size_t ind = 0; ind < huff_table_actual.size(); ind++)
        {
            REQUIRE( huff_table_actual[ind] == Bit_String{} );
        }
    }
    SECTION( "construction using initializer list" ) {
        using namespace std::string_literals;
        std::vector<Bit_String> expected_bit_strings{
            "0110"s,
            ""s,
            "11"s,
            "101"s,
            "0"s
        };

        Huff_Table huff_table_actual{ 
            "0110"s,
            ""s,
            "11"s,
            "101"s,
            "0"s, 
        };

        REQUIRE( huff_table_actual.size() == expected_size );

        // Check elements are default constructed
        for (size_t ind = 0; ind < huff_table_actual.size(); ind++)
        {
            REQUIRE( huff_table_actual[ind] == expected_bit_strings[ind] );
        }
    }
}

TEST_CASE( "Huff_Table::Load DC tables", "[Huff_Table]" ) {
    size_t expected_size{ 12 };

    SECTION( "Luminance tables" ) {
        Huff_Table huff_table_actual{ Huff_Table::load_DC_table(Image_Component::Luminance) };

        REQUIRE( huff_table_actual.size() == expected_size );

        // Check random selection of elements
        REQUIRE( huff_table_actual[0] == Bit_String{"00"} );
        REQUIRE( huff_table_actual[4] == Bit_String{"101"} );
        REQUIRE( huff_table_actual[11] == Bit_String{"111111110"} );
    }
    SECTION( "Chromiance tables" ) {
        Huff_Table huff_table_actual{ Huff_Table::load_DC_table(Image_Component::Chrominance) };

        REQUIRE( huff_table_actual.size() == expected_size );

        // Check random selection of elements
        REQUIRE( huff_table_actual[0] == Bit_String{"00"} );
        REQUIRE( huff_table_actual[4] == Bit_String{"1110"} );
        REQUIRE( huff_table_actual[11] == Bit_String{"11111111110"} );
    }
}

TEST_CASE( "Huff_Table::Load AC tables", "[Huff_Table]" ) {
    size_t expected_size{ 251 };

    SECTION( "Luminance tables" ) {
        Huff_Table huff_table_actual{ Huff_Table::load_AC_table(Image_Component::Luminance) };

        REQUIRE( huff_table_actual.size() == expected_size );

        // Check random selection of elements
        REQUIRE( huff_table_actual[0] == Bit_String{"1010"} );
        REQUIRE( huff_table_actual[0x32] == Bit_String{"111110111"} );
        REQUIRE( huff_table_actual[0x77] == Bit_String{"1111111110110010"} );
        REQUIRE( huff_table_actual[0xA1] == Bit_String{"111111010"} );
        REQUIRE( huff_table_actual[0xC5] == Bit_String{"1111111111011100"} );
        REQUIRE( huff_table_actual[0xFA] == Bit_String{"1111111111111110"} );
    }
    SECTION( "Chromiance tables" ) {
        Huff_Table huff_table_actual{ Huff_Table::load_AC_table(Image_Component::Chrominance) };

        REQUIRE( huff_table_actual.size() == expected_size );

        // Check random selection of elements
        REQUIRE( huff_table_actual[0] == Bit_String{"00"} );
        REQUIRE( huff_table_actual[0x32] == Bit_String{"11111000"} );
        REQUIRE( huff_table_actual[0x77] == Bit_String{"1111111110110011"} );
        REQUIRE( huff_table_actual[0xA1] == Bit_String{"111111000"} );
        REQUIRE( huff_table_actual[0xC5] == Bit_String{"1111111111011110"} );
        REQUIRE( huff_table_actual[0xFA] == Bit_String{"1111111111111110"} );
    }
}

