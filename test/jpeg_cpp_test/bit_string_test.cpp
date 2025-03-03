#include <catch2/catch_test_macros.hpp>
#include "jpeg_cpp/bit_string.h"

#include <sstream>
#include <string>

using namespace JPEG;

TEST_CASE( "Bit_String::byte iteration - integer bytes", "[Bit_String]" ) {
    Bit_String bs{ "1001010101100111" };
    std::vector<unsigned char> expected_result{
        0b10010101, 0b01100111
    };

    size_t ind{ 0 };

    for (auto byte_ptr=bs.begin_bytes(); byte_ptr!=bs.end_bytes(); ++byte_ptr)
    {
        // Double check we haven't gotten more bytes that we expected
        REQUIRE( ind < expected_result.size() );
        CHECK( *byte_ptr == expected_result[ind] );

        ind++;
    }

    // Double check we got all the bytes we expected
    CHECK( ind == expected_result.size() );
}

TEST_CASE( "Bit_String::byte iteration - non-intger bytes", "[Bit_String]" ) {
    Bit_String bs{ "1001010101100" };
    std::vector<unsigned char> expected_result{
        0b10010101, 0b01100000
    };

    size_t ind{ 0 };

    for (auto byte_ptr=bs.begin_bytes(); byte_ptr!=bs.end_bytes(); ++byte_ptr)
    {
        // Double check we haven't gotten more bytes that we expected
        REQUIRE( ind < expected_result.size() );
        CHECK( *byte_ptr == expected_result[ind] );

        ind++;
    }

    // Double check we got all the bytes we expected
    CHECK( ind == expected_result.size() );
}

TEST_CASE( "Bit_String::equality of Bit_Strings", "[Bit_String]" ) {
    Bit_String bs1{ 10 };

    // Set some bits
    bs1[0] = 1;
    bs1[4] = 1;
    bs1[9] = 1;

    REQUIRE( bs1.size() == 10 );

    SECTION( "ensure Bit_Strings are equal with themselves" ) {
        REQUIRE( bs1 == bs1 );
    }
    SECTION( "Bit_Strings of different lengths are not equal" ) {
        Bit_String bs2{ 9 };

        REQUIRE( bs1 != bs2 );
    }
    SECTION( "Bit_Strings of same length but with different bits are not equal" ) {
        Bit_String bs2{ 9 };

        bs1[0] = 1;
        bs1[9] = 1;

        REQUIRE( bs1 != bs2 );
    }
}

TEST_CASE( "Bit_String::bits can be appended", "[Bit_String]" ) {

    SECTION( "when Bit_String::size() is a multiple of 8" ) {
        Bit_String bs_initial{ "" };
        Bit_String bs_expected{ "1" };

        bs_initial.append_bit(1);

        REQUIRE( bs_initial == bs_expected );
    }
    SECTION( "when Bit_String::size() is not a multiple of 8" ) {
        Bit_String bs_initial{ "100" };
        Bit_String bs_expected{ "1001" };

        bs_initial.append_bit(1);

        REQUIRE( bs_initial == bs_expected );
    }
}

TEST_CASE( "Bit_String::last ssss bits of unsigned int can be appended", "[Bit_String]" ) {

    const unsigned int test_n{ 0b001111001010u };

    SECTION( "ssss == 0" ) {
        Bit_String bs_initial{ "" };
        Bit_String bs_expected{ "" };
        int ssss{ 0 };

        bs_initial.append_last_ssss_bits(test_n, ssss);

        REQUIRE( bs_initial == bs_expected );
    }
    SECTION( "0 <= ssss <= 7" ) {
        Bit_String bs_initial{ "" };
        Bit_String bs_expected{ "01010" };
        int ssss{ 5 };

        bs_initial.append_last_ssss_bits(test_n, ssss);

        REQUIRE( bs_initial == bs_expected );
    }
    SECTION( "ssss > 8" ) {
        Bit_String bs_initial{ "" };
        Bit_String bs_expected{ "111001010" };
        int ssss{ 9 };

        bs_initial.append_last_ssss_bits(test_n, ssss);

        REQUIRE( bs_initial == bs_expected );
    }
}

TEST_CASE( "Bit_String::bytes can be appended", "[Bit_String]" ) {

    const unsigned char test_byte{ 0b11100101u };

    SECTION( "Bit_string::size() is a multiple of 8" ) {
        Bit_String bs_initial{ "" };
        Bit_String bs_expected{ "11100101" };

        bs_initial.append_byte(test_byte);

        REQUIRE( bs_initial == bs_expected );
    }
    SECTION( "Bit_string::size() is not a multiple of 8" ) {
        Bit_String bs_initial{ "00" };
        Bit_String bs_expected{ "0011100101" };

        bs_initial.append_byte(test_byte);
        
        REQUIRE( bs_initial == bs_expected );
    }
}

TEST_CASE( "Bit_String::can be extended with other Bit_Strings", "[Bit_String]" ) {

    Bit_String bs_for_extending{ "11100101111" };

    SECTION( "Bit_string::size() is a multiple of 8" ) {
        Bit_String bs_initial{ "" };
        Bit_String bs_expected{ "11100101111" };

        bs_initial.extend(bs_for_extending);

        REQUIRE( bs_initial == bs_expected );
    }
    SECTION( "Bit_string::size() is not a multiple of 8" ) {
        Bit_String bs_initial{ "00" };
        Bit_String bs_expected{ "0011100101111" };

        bs_initial.extend(bs_for_extending);
        
        REQUIRE( bs_initial == bs_expected );
    }
}

TEST_CASE( "Bit_String::Bit_Strings can be created from strings", "[Bit_String]" ) {

    SECTION( "empty string produces Bit_string of zero length" ) {
        Bit_String bs{ "" };

        REQUIRE( bs.size() == 0 );
    }
    SECTION( "string of non-zero length" ) {
        Bit_String bs_from_string{ "01001" };
        Bit_String bs_expected{ 5 };

        bs_expected[1] = 1;
        bs_expected[4] = 1;

        REQUIRE( bs_from_string == bs_expected );
    }
}

TEST_CASE( "Bit_String::oputput to iostream", "[Bit_String]" ) {

    SECTION( "empty Bit_string produces empty string" ) {
        Bit_String bs{ "" };
        std::string str_expected{};

        std::stringstream str_stream;

        str_stream << bs;

        std::string str_actual{ str_stream.str() };
        
        REQUIRE( str_actual == str_expected );
    }
    SECTION( "non-zero length Bit_String produces correct output" ) {
        Bit_String bs{ "1011" };
        std::string str_expected{ "1011" };

        std::stringstream str_stream;

        str_stream << bs;

        std::string str_actual{ str_stream.str() };
        
        REQUIRE( str_actual == str_expected );
    }
}
