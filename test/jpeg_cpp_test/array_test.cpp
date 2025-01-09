#include <catch2/catch_test_macros.hpp>
#include "jpeg_cpp/array.h"

#include <sstream>
#include <string>
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

TEST_CASE( "DU_Array::construction of DU_Array", "[DU_Array]" ) {

    SECTION( "constructing using shape" ) {
        size_t expected_du{ 2 }, expected_rows{ 8 }, expected_cols{ 8 };
        DU_Array<int> arr_actual{ 2 };

        REQUIRE( arr_actual.size() == expected_du*expected_rows*expected_cols );
        REQUIRE( arr_actual.shape()[0] == expected_du );
        REQUIRE( arr_actual.shape()[1] == expected_rows );
        REQUIRE( arr_actual.shape()[2] == expected_cols );
    }
    SECTION( "constructing using initializer list" ) {
        size_t expected_du{ 2 }, expected_rows{ 8 }, expected_cols{ 8 };
        DU_Array<int> arr_actual = {
            {
                {  0,   1,   2,   3,   4,   5,   6,   7},
                {  8,   9,  10,  11,  12,  13,  14,  15},
                { 16,  17,  18,  19,  20,  21,  22,  23},
                { 24,  25,  26,  27,  28,  29,  30,  31},
                { 32,  33,  34,  35,  36,  37,  38,  39},
                { 40,  41,  42,  43,  44,  45,  46,  47}, 
                { 48,  49,  50,  51,  52,  53,  54,  55}, 
                { 56,  57,  58,  59,  60,  61,  62,  63},
            },
            {
                { 64,  65,  66,  67,  68,  69,  70,  71}, 
                { 72,  73,  74,  75,  76,  77,  78,  79}, 
                { 80,  81,  82,  83,  84,  85,  86,  87}, 
                { 88,  89,  90,  91,  92,  93,  94,  95}, 
                { 96,  97,  98,  99, 100, 101, 102, 103},
                {104, 105, 106, 107, 108, 109, 110, 111},
                {112, 113, 114, 115, 116, 117, 118, 119},
                {120, 121, 122, 123, 124, 125, 126, 127}
            }
        };

        REQUIRE( arr_actual.size() == expected_du*expected_rows*expected_cols );
        REQUIRE( arr_actual.shape()[0] == expected_du );
        REQUIRE( arr_actual.shape()[1] == expected_rows );
        REQUIRE( arr_actual.shape()[2] == expected_cols );

        // Check values of elements
        for (size_t du_ind = 0; du_ind < expected_du; du_ind++)
        {
            for (size_t i = 0; i < expected_rows; i++)
            {
                for (size_t j = 0; j < expected_cols; j++)
                {
                    REQUIRE( arr_actual.at(du_ind, i, j) == du_ind*expected_rows*expected_cols + i*expected_cols + j );
                }
            }
        }
    }
}

TEST_CASE( "DU_Array::3d subscripting of DU_Array", "[DU_Array]" ) {
    size_t expected_du{ 2 }, expected_rows{ 8 }, expected_cols{ 8 };
    DU_Array<int> arr_actual{ expected_du };

    for (int i = 0; i < arr_actual.size(); i++)
    {
        arr_actual[i] = i;
    }

    // Now test each by accessing them in row major order
    // Check values of elements
        for (size_t du_ind = 0; du_ind < expected_du; du_ind++)
        {
            for (size_t i = 0; i < expected_rows; i++)
            {
                for (size_t j = 0; j < expected_cols; j++)
                {
                    REQUIRE( arr_actual.at(du_ind, i, j) == du_ind*expected_rows*expected_cols + i*expected_cols + j );
                }
            }
        }
}

TEST_CASE( "DU_Array::output to iostream", "[DU_Array]" ) {
    DU_Array<int> arr_input = {
        {
            {  0,   1,   2,   3,   4,   5,   6,   7},
            {  8,   9,  10,  11,  12,  13,  14,  15},
            { 16,  17,  18,  19,  20,  21,  22,  23},
            { 24,  25,  26,  27,  28,  29,  30,  31},
            { 32,  33,  34,  35,  36,  37,  38,  39},
            { 40,  41,  42,  43,  44,  45,  46,  47}, 
            { 48,  49,  50,  51,  52,  53,  54,  55}, 
            { 56,  57,  58,  59,  60,  61,  62,  63},
        },
        {
            { 64,  65,  66,  67,  68,  69,  70,  71}, 
            { 72,  73,  74,  75,  76,  77,  78,  79}, 
            { 80,  81,  82,  83,  84,  85,  86,  87}, 
            { 88,  89,  90,  91,  92,  93,  94,  95}, 
            { 96,  97,  98,  99, 100, 101, 102, 103},
            {104, 105, 106, 107, 108, 109, 110, 111},
            {112, 113, 114, 115, 116, 117, 118, 119},
            {120, 121, 122, 123, 124, 125, 126, 127}
        }
    };

    // Now create a string representation and check whether it is what we expect
    std::string str_expected{
        "[\n"
        "  [\n"
        "    [ 0, 1, 2, 3, 4, 5, 6, 7 ],\n"
        "    [ 8, 9, 10, 11, 12, 13, 14, 15 ],\n"
        "    [ 16, 17, 18, 19, 20, 21, 22, 23 ],\n"
        "    [ 24, 25, 26, 27, 28, 29, 30, 31 ],\n"
        "    [ 32, 33, 34, 35, 36, 37, 38, 39 ],\n"
        "    [ 40, 41, 42, 43, 44, 45, 46, 47 ],\n"
        "    [ 48, 49, 50, 51, 52, 53, 54, 55 ],\n"
        "    [ 56, 57, 58, 59, 60, 61, 62, 63 ],\n"
        "  ],\n"
        "  [\n"
        "    [ 64, 65, 66, 67, 68, 69, 70, 71 ],\n"
        "    [ 72, 73, 74, 75, 76, 77, 78, 79 ],\n"
        "    [ 80, 81, 82, 83, 84, 85, 86, 87 ],\n"
        "    [ 88, 89, 90, 91, 92, 93, 94, 95 ],\n"
        "    [ 96, 97, 98, 99, 100, 101, 102, 103 ],\n"
        "    [ 104, 105, 106, 107, 108, 109, 110, 111 ],\n"
        "    [ 112, 113, 114, 115, 116, 117, 118, 119 ],\n"
        "    [ 120, 121, 122, 123, 124, 125, 126, 127 ],\n"
        "  ]\n"
        "]"
    };

    // Compute the actual string and see if their the same
    std::stringstream str_stream;

    str_stream << arr_input;

    std::string str_actual{ str_stream.str() };
    
    REQUIRE( str_actual == str_expected );
    
}

