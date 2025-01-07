#include <catch2/catch_test_macros.hpp>
#include "jpeg_cpp/array.h"

#include <vector>

using namespace JPEG;

TEST_CASE( "Array::construction of Arrays", "[Array]" ) {

    SECTION( "can produce array of set size" ) {
        size_t expected_size{ 5 };
        Array<int> arr_actual(5);

        REQUIRE( arr_actual.size() == expected_size );
    }
    SECTION( "can produce array of set size with default constructed elements" ) {
        size_t expected_size{ 5 };
        std::vector<int> arr_expected = {1, 2, 3, 4, 5};
        Array<int> arr_actual = {1, 2, 3, 4, 5};

        REQUIRE( arr_actual.size() == expected_size );

        for (size_t i = 0; i < expected_size; i++)
        {
            REQUIRE( arr_actual[i] == arr_expected[i] );
        }
    }
}
