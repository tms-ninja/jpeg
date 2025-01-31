#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include <array>

#include "jpeg_cpp/array.h"
#include "jpeg_cpp/jpeg.h"

using namespace JPEG;
using namespace Catch::Matchers;


TEST_CASE( "encode_greyscale_image()::markers appear in correct order", "[encode_greyscale_image()]" ) {
    Array_2d<double> input_array{ 
        { 0, 1, 2, 3 },
        { 5, 6, 7, 8 }
     };

    const std::vector<unsigned char> encoded_image{ encode_greyscale_image(input_array) };

    // Now check order of markers
    // Must start with SOI, end with EOI
    // SOI and EOI must not appear anywhere else
    // Contain at least 1 DQT, 1 DHT, 1 SOF0, 1 SOS
    // SOS must come after SOF0
    // First check the encoded image is long enough to have all 6 markers

    REQUIRE( encoded_image.size() >= 12 );

    // Check for SOI at start
    REQUIRE( encoded_image[0] == 0xFF );
    REQUIRE( encoded_image[1] == 0xD8 );

    // Check for EOI at end
    REQUIRE( encoded_image[encoded_image.size()-2] == 0xFF );
    REQUIRE( encoded_image[encoded_image.size()-1] == 0xD9 );

    // Now check each marker in turn
    bool found_marker{ false };
    bool found_SOF{ false };

    std::array<unsigned char, 4> expected_markers{
        0xDB,   // DQT
        0xC4,   // DHT
        0xC0,   // SOF0
        0xDA,   // SOS
    };

    for (size_t ind = 2; ind < encoded_image.size()-2; ind++)
    {
        unsigned char cur_byte{ encoded_image[ind] };

        if (cur_byte==0xFF)
        {
            found_marker = true;
        }
        else if (found_marker)
        {
            found_marker = false;

            // Ignore byte stuffing
            if (cur_byte==0x00)
            {
                continue;
            }

            // Ensure it's one of the expected markers
            CAPTURE( ind );
            CHECK_THAT( expected_markers, Contains(cur_byte) );

            // Check for SOF0 and SOS
            if (cur_byte==0xC0)
            {
                found_SOF = true;
            }
            else if (cur_byte==0xDA)
            {
                // SOS marker, must follow a SOF marker
                CHECK( found_SOF );
            }
        }
    }
}
