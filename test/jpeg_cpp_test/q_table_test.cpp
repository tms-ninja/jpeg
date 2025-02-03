#include <catch2/catch_test_macros.hpp>
#include "jpeg_cpp/q_table.h"

#include <array>
#include <sstream>
#include <string>

using namespace JPEG;

TEST_CASE( "Q_Table::construction of Q_Tables", "[Q_Table]" ) {
    size_t expected_size{ 64 };
    std::array<size_t, 2> expected_shape{ 8, 8 };

    SECTION( "default construction" ) {
        Q_Table q_table_actual;

        REQUIRE( q_table_actual.size() == expected_size );
        REQUIRE( q_table_actual.shape()[0] == expected_shape[0] );
        REQUIRE( q_table_actual.shape()[1] == expected_shape[1] );
    }
    SECTION( "construction from initializer list" ) {
        Q_Table q_table_actual = {
            { 0,  1,  2,  3,  4,  5,  6,  7},
            { 8,  9, 10, 11, 12, 13, 14, 15}, 
            {16, 17, 18, 19, 20, 21, 22, 23},
            {24, 25, 26, 27, 28, 29, 30, 31},
            {32, 33, 34, 35, 36, 37, 38, 39},
            {40, 41, 42, 43, 44, 45, 46, 47},
            {48, 49, 50, 51, 52, 53, 54, 55},
            {56, 57, 58, 59, 60, 61, 62, 63}
        };

        REQUIRE( q_table_actual.size() == expected_size );
        REQUIRE( q_table_actual.shape()[0] == expected_shape[0] );
        REQUIRE( q_table_actual.shape()[1] == expected_shape[1] );

        for (size_t ind = 0; ind < expected_size; ind++)
        {
            CAPTURE( ind );
            CHECK( q_table_actual[ind] == ind );
        }
    }
}

TEST_CASE( "Q_Table::string output", "[Q_Table]" ) {
    Q_Table q_table = {
        { 0,  1,  2,  3,  4,  5,  6,  7},
        { 8,  9, 10, 11, 12, 13, 14, 15}, 
        {16, 17, 18, 19, 20, 21, 22, 23},
        {24, 25, 26, 27, 28, 29, 30, 31},
        {32, 33, 34, 35, 36, 37, 38, 39},
        {40, 41, 42, 43, 44, 45, 46, 47},
        {48, 49, 50, 51, 52, 53, 54, 55},
        {56, 57, 58, 59, 60, 61, 62, 63}
    };

    std::string str_expected{
        "[\n"
        "  [ 0, 1, 2, 3, 4, 5, 6, 7 ],\n"
        "  [ 8, 9, 10, 11, 12, 13, 14, 15 ],\n"
        "  [ 16, 17, 18, 19, 20, 21, 22, 23 ],\n"
        "  [ 24, 25, 26, 27, 28, 29, 30, 31 ],\n"
        "  [ 32, 33, 34, 35, 36, 37, 38, 39 ],\n"
        "  [ 40, 41, 42, 43, 44, 45, 46, 47 ],\n"
        "  [ 48, 49, 50, 51, 52, 53, 54, 55 ],\n"
        "  [ 56, 57, 58, 59, 60, 61, 62, 63 ]\n"
        "]"
    };

    // Compute the actual string and see if their the same
    std::stringstream str_stream;

    str_stream << q_table;

    std::string str_actual{ str_stream.str() };
    
    REQUIRE( str_actual == str_expected );
}

TEST_CASE( "Q_Table::generation from quality factor", "[Q_Table]" ) {
    Image_Component comp{ Image_Component::Luminance };
    
    SECTION("qf = 0") {
        int qf{ 0 };
        Q_Table table{ Q_Table::load_q_table_from_quality_factor(comp, qf) };

        CAPTURE( qf );

        // quantization table elements of 8 bit precision images
        // can't be more than 8 bits
        CHECK( table(0, 0) == 255 );
        CHECK( table(4, 0) == 255 );
        CHECK( table(0, 4) == 255 );
        CHECK( table(6, 7) == 255 );
    }
    SECTION("qf = 25") {
        int qf{ 25 };
        Q_Table table{ Q_Table::load_q_table_from_quality_factor(comp, qf) };

        CAPTURE( qf );

        CHECK( table(0, 0) == 32 );
        CHECK( table(4, 0) == 36 );
        CHECK( table(0, 4) == 48 );
        CHECK( table(6, 7) == 202 );
    }
    SECTION("qf = 50") {
        // Should generate spec tables
        int qf{ 50 };
        Q_Table table{ Q_Table::load_q_table_from_quality_factor(comp, qf) };

        CAPTURE( qf );

        CHECK( table(0, 0) == 16 );
        CHECK( table(4, 0) == 18 );
        CHECK( table(0, 4) == 24 );
        CHECK( table(6, 7) == 101 );
    }
    SECTION("qf = 95") {
        // Should generate spec tables
        int qf{ 95 };
        Q_Table table{ Q_Table::load_q_table_from_quality_factor(comp, qf) };

        CAPTURE( qf );

        CHECK( table(0, 0) == 2 );
        CHECK( table(4, 0) == 2 );
        CHECK( table(0, 4) == 2 );
        CHECK( table(6, 7) == 10 );
    }
    SECTION("qf = 100") {
        // Should generate spec tables
        int qf{ 100 };
        Q_Table table{ Q_Table::load_q_table_from_quality_factor(comp, qf) };

        CAPTURE( qf );

        CHECK( table(0, 0) == 1 );
        CHECK( table(4, 0) == 1 );
        CHECK( table(0, 4) == 1 );
        CHECK( table(6, 7) == 1 );
    }
    SECTION("qf < 0 rejected") {
        // Should generate spec tables
        int qf{ -1 };

        CHECK_THROWS_AS( Q_Table::load_q_table_from_quality_factor(comp, qf), std::invalid_argument );
    }
    SECTION("qf > 100 rejected") {
        // Should generate spec tables
        int qf{ 101 };

        CHECK_THROWS_AS( Q_Table::load_q_table_from_quality_factor(comp, qf), std::invalid_argument );
    }
}
