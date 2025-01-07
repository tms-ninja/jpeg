#include <catch2/catch_test_macros.hpp>
#include "jpeg_cpp/bit_string.h"


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
        Bit_String bs_initial{ Bit_String::from_string("") };
        Bit_String bs_expected{ Bit_String::from_string("1") };

        bs_initial.append_bit(1);

        REQUIRE( bs_initial == bs_expected );
    }
    SECTION( "when Bit_String::size() is not a multiple of 8" ) {
        Bit_String bs_initial{ Bit_String::from_string("100") };
        Bit_String bs_expected{ Bit_String::from_string("1001") };

        bs_initial.append_bit(1);

        REQUIRE( bs_initial == bs_expected );
    }
}

TEST_CASE( "Bit_String::last ssss bits of unsigned int can be appended", "[Bit_String]" ) {

    const unsigned int test_n{ 0b001111001010u };

    SECTION( "ssss == 0" ) {
        Bit_String bs_initial{ Bit_String::from_string("") };
        Bit_String bs_expected{ Bit_String::from_string("") };
        int ssss{ 0 };

        bs_initial.append_last_ssss_bits(test_n, ssss);

        REQUIRE( bs_initial == bs_expected );
    }
    SECTION( "0 <= ssss <= 7" ) {
        Bit_String bs_initial{ Bit_String::from_string("") };
        Bit_String bs_expected{ Bit_String::from_string("01010") };
        int ssss{ 5 };

        bs_initial.append_last_ssss_bits(test_n, ssss);

        REQUIRE( bs_initial == bs_expected );
    }
    SECTION( "ssss > 8" ) {
        Bit_String bs_initial{ Bit_String::from_string("") };
        Bit_String bs_expected{ Bit_String::from_string("111001010") };
        int ssss{ 9 };

        bs_initial.append_last_ssss_bits(test_n, ssss);

        REQUIRE( bs_initial == bs_expected );
    }
}

TEST_CASE( "Bit_String::bytes can be appended", "[Bit_String]" ) {

    const unsigned char test_byte{ 0b11100101u };

    SECTION( "Bit_string::size() is a multiple of 8" ) {
        Bit_String bs_initial{ Bit_String::from_string("") };
        Bit_String bs_expected{ Bit_String::from_string("11100101") };

        bs_initial.append_byte(test_byte);

        REQUIRE( bs_initial == bs_expected );
    }
    SECTION( "Bit_string::size() is not a multiple of 8" ) {
        Bit_String bs_initial{ Bit_String::from_string("00") };
        Bit_String bs_expected{ Bit_String::from_string("0011100101") };

        bs_initial.append_byte(test_byte);
        
        REQUIRE( bs_initial == bs_expected );
    }
}

TEST_CASE( "Bit_String::can be extended with other Bit_Strings", "[Bit_String]" ) {

    Bit_String bs_for_extending{ Bit_String::from_string("11100101111") };

    SECTION( "Bit_string::size() is a multiple of 8" ) {
        Bit_String bs_initial{ Bit_String::from_string("") };
        Bit_String bs_expected{ Bit_String::from_string("11100101111") };

        bs_initial.extend(bs_for_extending);

        REQUIRE( bs_initial == bs_expected );
    }
    SECTION( "Bit_string::size() is not a multiple of 8" ) {
        Bit_String bs_initial{ Bit_String::from_string("00") };
        Bit_String bs_expected{ Bit_String::from_string("0011100101111") };

        bs_initial.extend(bs_for_extending);
        
        REQUIRE( bs_initial == bs_expected );
    }
}

TEST_CASE( "Bit_String::Bit_Strings can be created from strings", "[Bit_String]" ) {

    SECTION( "empty string produces Bit_string of zero length" ) {
        Bit_String bs{ Bit_String::from_string("") };

        REQUIRE( bs.size() == 0 );
    }
    SECTION( "string of non-zero length" ) {
        Bit_String bs_from_string{ Bit_String::from_string("01001") };
        Bit_String bs_expected{ 5 };

        bs_expected[1] = 1;
        bs_expected[4] = 1;

        REQUIRE( bs_from_string == bs_expected );
    }
}


