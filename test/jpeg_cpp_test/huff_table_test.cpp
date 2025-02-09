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
            CAPTURE( ind );
            CHECK( huff_table_actual[ind] == Bit_String{} );
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
            CAPTURE( ind );
            CHECK( huff_table_actual[ind] == expected_bit_strings[ind] );
        }
    }
}

TEST_CASE( "Huff_Table::Load DC tables", "[Huff_Table]" ) {
    size_t expected_size{ 12 };

    SECTION( "Luminance tables" ) {
        Huff_Table huff_table_actual{ Huff_Table::load_DC_table(Image_Component::Luminance) };

        REQUIRE( huff_table_actual.size() == expected_size );

        // Check random selection of elements
        CHECK( huff_table_actual[0] == Bit_String{"00"} );
        CHECK( huff_table_actual[4] == Bit_String{"101"} );
        CHECK( huff_table_actual[11] == Bit_String{"111111110"} );
    }
    SECTION( "Chromiance tables" ) {
        Huff_Table huff_table_actual{ Huff_Table::load_DC_table(Image_Component::Chrominance) };

        REQUIRE( huff_table_actual.size() == expected_size );

        // Check random selection of elements
        CHECK( huff_table_actual[0] == Bit_String{"00"} );
        CHECK( huff_table_actual[4] == Bit_String{"1110"} );
        CHECK( huff_table_actual[11] == Bit_String{"11111111110"} );
    }
}

TEST_CASE( "Huff_Table::Load AC tables", "[Huff_Table]" ) {
    size_t expected_size{ 251 };

    SECTION( "Luminance tables" ) {
        Huff_Table huff_table_actual{ Huff_Table::load_AC_table(Image_Component::Luminance) };

        REQUIRE( huff_table_actual.size() == expected_size );

        // Check random selection of elements
        CHECK( huff_table_actual[0] == Bit_String{"1010"} );
        CHECK( huff_table_actual[0x32] == Bit_String{"111110111"} );
        CHECK( huff_table_actual[0x77] == Bit_String{"1111111110110010"} );
        CHECK( huff_table_actual[0xA1] == Bit_String{"111111010"} );
        CHECK( huff_table_actual[0xC5] == Bit_String{"1111111111011100"} );
        CHECK( huff_table_actual[0xFA] == Bit_String{"1111111111111110"} );
    }
    SECTION( "Chromiance tables" ) {
        Huff_Table huff_table_actual{ Huff_Table::load_AC_table(Image_Component::Chrominance) };

        REQUIRE( huff_table_actual.size() == expected_size );

        // Check random selection of elements
        CHECK( huff_table_actual[0] == Bit_String{"00"} );
        CHECK( huff_table_actual[0x32] == Bit_String{"11111000"} );
        CHECK( huff_table_actual[0x77] == Bit_String{"1111111110110011"} );
        CHECK( huff_table_actual[0xA1] == Bit_String{"111111000"} );
        CHECK( huff_table_actual[0xC5] == Bit_String{"1111111111011110"} );
        CHECK( huff_table_actual[0xFA] == Bit_String{"1111111111111110"} );
    }
}

TEST_CASE( "Huff_Table::load_table_from_BITS_and_HUFFVAL()", "[Huff_Table]" ) {
    std::array<unsigned int, 16> bits{
        0, 1, 3, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };
    std::vector<unsigned int> huffval{
        3, 6, 1, 0, 9, 7
    };
    std::vector<Bit_String> expected_codes{
        Bit_String("00"),
        Bit_String("010"),
        Bit_String("011"),
        Bit_String("100"),
        Bit_String("10100"),
        Bit_String("10101")
    };
    
    Huff_Table actual_result{ Huff_Table::load_table_from_BITS_and_HUFFVAL(bits, huffval) };

    // Maximum symbol in HUFFVAL is 9, so there must be at least 10 entries in the 
    // Huffman table. (Note that not all of these entries necessarily have a Huffman 
    // code)
    REQUIRE( actual_result.size()>= 10 );

    // Check each symbol has the expected code
    for (size_t symbol_ind = 0; symbol_ind < huffval.size(); symbol_ind++)
    {
        CAPTURE( symbol_ind );
        CHECK( actual_result[huffval[symbol_ind]] == expected_codes[symbol_ind] );
    }
}

TEST_CASE( "Huff_Table::gen_table_from_stats()", "[Huff_Table]" ) {
    std::vector<unsigned int> stats{
        4,
        7,
        4,
        0,
        3
    };

    Huff_Table expected_result(256);

    expected_result[0] = Bit_String{"00"};
    expected_result[1] = Bit_String{"01"};
    expected_result[2] = Bit_String{"10"};
    expected_result[3] = Bit_String{};
    expected_result[4] = Bit_String{"110"};

    Huff_Table actual_result{ Huff_Table::gen_table_from_stats(stats) };

    REQUIRE( actual_result.size() == expected_result.size() );

    // Check each code one by one
    for (size_t code_ind = 0; code_ind < actual_result.size(); code_ind++)
    {
        CAPTURE( code_ind );
        CHECK( actual_result[code_ind] == expected_result[code_ind] );
    }
}