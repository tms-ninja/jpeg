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

TEST_CASE( "Array::resizing of Arrays", "[Array]" ) {
    Array<int> arr_actual = {0, 1, 2, 3, 4, 5};
    size_t initial_size{ 5 };

    SECTION( "can be resized larger" ) {
        size_t new_size{ 10 };

        arr_actual.resize(new_size);

        // Check the new size
        REQUIRE( arr_actual.size() == new_size );

        // Check initial elements are preserved
        for (size_t ind = 0; ind < initial_size; ind++)
        {
            REQUIRE( arr_actual[ind] == ind );
        }
    }
    SECTION( "can be resized smaller" ) {
        size_t new_size{ 3 };

        arr_actual.resize(new_size);

        // Check the new size
        REQUIRE( arr_actual.size() == new_size );

        // Check initial elements are preserved
        for (size_t ind = 0; ind < new_size; ind++)
        {
            REQUIRE( arr_actual[ind] == ind );
        }
    }
    SECTION( "can be resized the same" ) {
        size_t new_size{ initial_size };

        arr_actual.resize(new_size);

        // Check the new size
        REQUIRE( arr_actual.size() == new_size );

        // Check initial elements are preserved
        for (size_t ind = 0; ind < new_size; ind++)
        {
            REQUIRE( arr_actual[ind] == ind );
        }
    }
}

TEST_CASE( "Array_2d::construction of Array_2ds", "[Array_2d]" ) {

    SECTION( "constructing using number of rows and columns" ) {
        size_t expected_rows{ 2 }, expected_cols{ 3 };
        Array_2d<int> arr_actual{ 2, 3 };

        REQUIRE( arr_actual.size() == expected_rows*expected_cols );
        REQUIRE( arr_actual.shape()[0] == expected_rows );
        REQUIRE( arr_actual.shape()[1] == expected_cols );
    }
    SECTION( "constructing using initializer list" ) {
        size_t expected_rows{ 2 }, expected_cols{ 3 };
        Array_2d<int> arr_actual = {
            {0, 1, 2},
            {3, 4, 5}
        };

        REQUIRE( arr_actual.size() == expected_rows*expected_cols );
        REQUIRE( arr_actual.shape()[0] == expected_rows );
        REQUIRE( arr_actual.shape()[1] == expected_cols );

        // Check values of elements
        for (size_t i = 0; i < expected_rows; i++)
        {
            for (size_t j = 0; j < expected_cols; j++)
            {
                REQUIRE( arr_actual.at(i, j) == i*expected_cols + j );
            }
        }
    }
}

TEST_CASE( "Array_2d::2d subscripting of Array_2ds", "[Array_2d]" ) {
    size_t expected_rows{ 2 }, expected_cols{ 3 };
    Array_2d<int> arr_actual{ 2, 3 };

    for (int i = 0; i < arr_actual.size(); i++)
    {
        arr_actual[i] = i;
    }

    // Now test each by accessing them in row major order
    for (size_t i = 0; i < expected_rows; i++)
    {
        for (size_t j = 0; j < expected_cols; j++)
        {
            REQUIRE( arr_actual.at(i, j) == i*expected_cols + j );
        }
    }
}

