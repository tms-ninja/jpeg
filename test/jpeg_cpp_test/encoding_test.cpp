#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
#include <stdexcept>

#include "jpeg_cpp/array.h"
#include "jpeg_cpp/encoding.h"
#include "jpeg_cpp/huff_table.h"

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

// Generates a DU_Array with a single data unit containing example on the 
// JPEG Wikipedia page
DU_Array<double> gen_wiki_DU_Array()
{
    DU_Array<double> array{{
        {52, 55, 61,  66,  70,  61, 64, 73},
        {63, 59, 55,  90, 109,  85, 69, 72},
        {62, 59, 68, 113, 144, 104, 66, 73},
        {63, 58, 71, 122, 154, 106, 70, 69},
        {67, 61, 68, 104, 126,  88, 68, 70},
        {79, 65, 60,  70,  77,  68, 58, 75},
        {85, 71, 64,  59,  55,  61, 65, 83},
        {87, 79, 69,  68,  65,  76, 78, 94}
    }};

    return array;
}

// Generates a DU_Array with two data units. The first is the Wiki example, the second
// its transpose
DU_Array<double> gen_2_DU_Array()
{
    DU_Array<double> array{
        {
            {52, 55, 61,  66,  70,  61, 64, 73},
            {63, 59, 55,  90, 109,  85, 69, 72},
            {62, 59, 68, 113, 144, 104, 66, 73},
            {63, 58, 71, 122, 154, 106, 70, 69},
            {67, 61, 68, 104, 126,  88, 68, 70},
            {79, 65, 60,  70,  77,  68, 58, 75},
            {85, 71, 64,  59,  55,  61, 65, 83},
            {87, 79, 69,  68,  65,  76, 78, 94}
        },
        {
            {52,  63,  62,  63,  67,  79,  85,  87},
            {55,  59,  59,  58,  61,  65,  71,  79},
            {61,  55,  68,  71,  68,  60,  64,  69},
            {66,  90, 113, 122, 104,  70,  59,  68},
            {70, 109, 144, 154, 126,  77,  55,  65},
            {61,  85, 104, 106,  88,  68,  61,  76},
            {64,  69,  66,  70,  68,  58,  65,  78},
            {73,  72,  73,  69,  70,  75,  83,  94}
        }
    };

    return array;
}

// Loads the quantization table on the Wikipedia JPEG page
Q_Table gen_wiki_q_table()
{
    Q_Table table{
        {16, 11, 10, 16, 24, 40, 51, 61},
        {12, 12, 14, 19, 26, 58, 60, 55},
        {14, 13, 16, 24, 40, 57, 69, 56},
        {14, 17, 22, 29, 51, 87, 80, 62},
        {18, 22, 37, 56, 68, 109, 103, 77},
        {24, 35, 55, 64, 81, 104, 113, 92},
        {49, 64, 78, 87, 103, 121, 120, 101},
        {72, 92, 95, 98, 112, 100, 103, 99}
    };

    return table;
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

TEST_CASE( "load_DCT_matrix()::loads the correct matrix", "[load_DCT_matrix()]" ) {
    Array_2d<double> array{ load_DCT_matrix() };

    REQUIRE( array.shape()[0] == 8 );
    REQUIRE( array.shape()[1] == 8 );
    
    // Check a selection of values
    REQUIRE_THAT( array(0, 0), WithinRel(0.5) );
    REQUIRE_THAT( array(0, 1), WithinRel(0.4903926402016152) );
    REQUIRE_THAT( array(1, 0), WithinRel(0.5) );
    REQUIRE_THAT( array(3, 4), WithinRel(0.3535533905932737) );
    REQUIRE_THAT( array(7, 2), WithinRel(0.4619397662556433) );
}

TEST_CASE( "load_DCT_matrix_transpose()::loads the correct matrix", "[load_DCT_matrix_transpose()]" ) {
    Array_2d<double> array{ load_DCT_matrix_transpose() };

    REQUIRE( array.shape()[0] == 8 );
    REQUIRE( array.shape()[1] == 8 );
    
    // Check a selection of values
    REQUIRE_THAT( array(0, 0), WithinRel(0.5) );
    REQUIRE_THAT( array(0, 1), WithinRel(0.5) );
    REQUIRE_THAT( array(1, 0), WithinRel(0.4903926402016152) );
    REQUIRE_THAT( array(3, 4), WithinRel(0.2777851165098009) );
    REQUIRE_THAT( array(7, 2), WithinRel(0.4157348061512727) );
}

TEST_CASE( "mat_mul()::multiplies matrices correctly", "[mat_mul()]" ) {
    Array_2d<double> mat_a{ 
        { 0,  1,  2,  3,  4,  5,  6,  7},
        { 8,  9, 10, 11, 12, 13, 14, 15},
        {16, 17, 18, 19, 20, 21, 22, 23},
        {24, 25, 26, 27, 28, 29, 30, 31},
        {32, 33, 34, 35, 36, 37, 38, 39},
        {40, 41, 42, 43, 44, 45, 46, 47},
        {48, 49, 50, 51, 52, 53, 54, 55},
        {56, 57, 58, 59, 60, 61, 62, 63}
    };

    Array_2d<double> mat_b{ 
        { 64,  65,  66,  67,  68,  69,  70,  71},
        { 72,  73,  74,  75,  76,  77,  78,  79},
        { 80,  81,  82,  83,  84,  85,  86,  87},
        { 88,  89,  90,  91,  92,  93,  94,  95},
        { 96,  97,  98,  99, 100, 101, 102, 103},
        {104, 105, 106, 107, 108, 109, 110, 111},
        {112, 113, 114, 115, 116, 117, 118, 119},
        {120, 121, 122, 123, 124, 125, 126, 127}
    };

    Array_2d<double> expected_result{ 
        { 2912,  2940,  2968,  2996,  3024,  3052,  3080,  3108},
        { 8800,  8892,  8984,  9076,  9168,  9260,  9352,  9444},
        {14688, 14844, 15000, 15156, 15312, 15468, 15624, 15780},
        {20576, 20796, 21016, 21236, 21456, 21676, 21896, 22116},
        {26464, 26748, 27032, 27316, 27600, 27884, 28168, 28452},
        {32352, 32700, 33048, 33396, 33744, 34092, 34440, 34788},
        {38240, 38652, 39064, 39476, 39888, 40300, 40712, 41124},
        {44128, 44604, 45080, 45556, 46032, 46508, 46984, 47460}
    };

    // Simple buffer to store the result
    std::array<double, 64> actual_result;

    mat_mul(mat_a.data(), mat_b.data(), &actual_result[0], mat_a.shape()[0]);

    for (size_t ind = 0; ind < actual_result.size(); ind++)
    {
        REQUIRE_THAT( actual_result[ind], WithinRel(expected_result[ind]) );
    }
}

TEST_CASE( "apply_DCT()::applies correct transform to a single data unit", "[apply_DCT()]" ) {
    DU_Array<double> input_array{ gen_wiki_DU_Array() };
    // Need to apply the level shift so we're consistent with Wikipedia
    apply_level_shift(input_array);

    // Note that as the factors of 1/sqrt(2) are absorbed into the quantization
    // step, the expected output is slightly different from what is on Wikipedia
    DU_Array<double> expected_result{{
        {-8.3075000000000000e+02, -4.2689050762823221e+01, -8.6545714987494691e+01, 3.8522219303704155e+01, 7.9372736188189933e+01, -2.8418867287080044e+01, -3.3766429043248785e+00, 6.5310566202395337e-01 },
        { 6.3152041816377666e+00, -2.1857439332259819e+01, -6.0758038116534024e+01, 1.0253636818417853e+01, 1.3145110120476241e+01, -7.0874180078452067e+00, -8.5354367129694957e+00, 4.8768884966804009e+00 },
        {-6.6233963509334060e+01, 7.3705973534266889e+00, 7.7129387578755527e+01, -2.4561982249733369e+01, -2.8911688429320645e+01, 9.9335209527750816e+00, 5.4168154723945428e+00, -5.6489508621374700e+00 },
        {-6.8638808107720052e+01, 1.2068360940019206e+01, 3.4099767172715055e+01, -1.4759411080801923e+01, -1.0240606801750433e+01, 6.2959674383730153e+00, 1.8311650530957329e+00, 1.9459365148648098e+00 },
        { 1.7147339443773728e+01, -6.5534499288920696e+00, -1.3196120970971871e+01, -3.9514277279078320e+00, -1.8749999999999907e+00, 1.7452844510267289e+00, -2.7872282503369412e+00, 3.1352823039767745e+00 },
        {-1.0938579410340271e+01, 2.9054613828905698e+00, 2.3797957648755750e+00, -5.9393139358655374e+00, -2.3777967067325991e+00, 9.4139159614138390e-01, 4.3037133436227535e+00, 1.8486910259091240e+00 },
        {-1.4575931682733301e+00, 1.8306744355203072e-01, 4.1681547239454325e-01, -2.4155613745353905e+00, -8.7779391994230582e-01, -3.0193065522845322e+00, 4.1206124212444859e+00, -6.6194845393858159e-01 },
        {-2.3387641928580027e-01, 1.4160712244183549e-01, -1.0715363895103391e+00, -4.1929120780447091e+00, -1.1703140920062558e+00, -9.7761079337536458e-02, 5.0126939164458373e-01, 1.6754588169203766e+00}
    }};

    apply_DCT(input_array);

    // Verify the DCT was performed correctly
    for (size_t ind = 0; ind < input_array.size(); ind++)
    {
        // Note sure why some seem so off, increase the tolerance for now
        REQUIRE_THAT( input_array[ind], WithinRel(expected_result[ind], 1e-12) );
    }
}

TEST_CASE( "apply_DCT()::applies correct transform to multiple data units", "[apply_DCT()]" ) {
    DU_Array<double> input_array{ gen_2_DU_Array() };
    // Need to apply the level shift so we're consistent with Wikipedia
    apply_level_shift(input_array);

    // Note that as the factors of 1/sqrt(2) are absorbed into the quantization
    // step, the expected output is slightly different from what is on Wikipedia
    DU_Array<double> expected_result{
        {
            {-8.3075000000000000e+02, -4.2689050762823221e+01, -8.6545714987494691e+01, 3.8522219303704155e+01, 7.9372736188189933e+01, -2.8418867287080044e+01, -3.3766429043248785e+00, 6.5310566202395337e-01 },
            { 6.3152041816377666e+00, -2.1857439332259819e+01, -6.0758038116534024e+01, 1.0253636818417853e+01, 1.3145110120476241e+01, -7.0874180078452067e+00, -8.5354367129694957e+00, 4.8768884966804009e+00 },
            {-6.6233963509334060e+01, 7.3705973534266889e+00, 7.7129387578755527e+01, -2.4561982249733369e+01, -2.8911688429320645e+01, 9.9335209527750816e+00, 5.4168154723945428e+00, -5.6489508621374700e+00 },
            {-6.8638808107720052e+01, 1.2068360940019206e+01, 3.4099767172715055e+01, -1.4759411080801923e+01, -1.0240606801750433e+01, 6.2959674383730153e+00, 1.8311650530957329e+00, 1.9459365148648098e+00 },
            { 1.7147339443773728e+01, -6.5534499288920696e+00, -1.3196120970971871e+01, -3.9514277279078320e+00, -1.8749999999999907e+00, 1.7452844510267289e+00, -2.7872282503369412e+00, 3.1352823039767745e+00 },
            {-1.0938579410340271e+01, 2.9054613828905698e+00, 2.3797957648755750e+00, -5.9393139358655374e+00, -2.3777967067325991e+00, 9.4139159614138390e-01, 4.3037133436227535e+00, 1.8486910259091240e+00 },
            {-1.4575931682733301e+00, 1.8306744355203072e-01, 4.1681547239454325e-01, -2.4155613745353905e+00, -8.7779391994230582e-01, -3.0193065522845322e+00, 4.1206124212444859e+00, -6.6194845393858159e-01 },
            {-2.3387641928580027e-01, 1.4160712244183549e-01, -1.0715363895103391e+00, -4.1929120780447091e+00, -1.1703140920062558e+00, -9.7761079337536458e-02, 5.0126939164458373e-01, 1.6754588169203766e+00}
        },
        {
            {-8.3075000000000000e+02, 6.3152041816377666e+00, -6.6233963509334060e+01, -6.8638808107720052e+01, 1.7147339443773728e+01, -1.0938579410340271e+01, -1.4575931682733301e+00, -2.3387641928580027e-01},
            {-4.2689050762823221e+01, -2.1857439332259819e+01, 7.3705973534266889e+00, 1.2068360940019206e+01, -6.5534499288920696e+00, 2.9054613828905698e+00, 1.8306744355203072e-01, 1.4160712244183549e-01},
            {-8.6545714987494691e+01, -6.0758038116534024e+01, 7.7129387578755527e+01, 3.4099767172715055e+01, -1.3196120970971871e+01, 2.3797957648755750e+00, 4.1681547239454325e-01, -1.0715363895103391e+00},
            { 3.8522219303704155e+01, 1.0253636818417853e+01, -2.4561982249733369e+01, -1.4759411080801923e+01, -3.9514277279078320e+00, -5.9393139358655374e+00, -2.4155613745353905e+00, -4.1929120780447091e+00},
            { 7.9372736188189933e+01, 1.3145110120476241e+01, -2.8911688429320645e+01, -1.0240606801750433e+01, -1.8749999999999907e+00, -2.3777967067325991e+00, -8.7779391994230582e-01, -1.1703140920062558e+00},
            {-2.8418867287080044e+01, -7.0874180078452067e+00, 9.9335209527750816e+00, 6.2959674383730153e+00, 1.7452844510267289e+00, 9.4139159614138390e-01, -3.0193065522845322e+00, -9.7761079337536458e-02},
            {-3.3766429043248785e+00, -8.5354367129694957e+00, 5.4168154723945428e+00, 1.8311650530957329e+00, -2.7872282503369412e+00, 4.3037133436227535e+00, 4.1206124212444859e+00, 5.0126939164458373e-01},
            { 6.5310566202395337e-01, 4.8768884966804009e+00, -5.6489508621374700e+00, 1.9459365148648098e+00, 3.1352823039767745e+00, 1.8486910259091240e+00, -6.6194845393858159e-01, 1.6754588169203766e+00}
        }
    };

    apply_DCT(input_array);

    // Verify the DCT was performed correctly
    for (size_t ind = 0; ind < input_array.size(); ind++)
    {
        // Note sure why some seem so off, increase the tolerance for now
        REQUIRE_THAT( input_array[ind], WithinRel(expected_result[ind], 1e-12) );
    }
}

TEST_CASE( "apply_quantization()::correctly quantizes a single data unit", "[apply_quantization()]" ) {
    DU_Array<double> input_array{ gen_wiki_DU_Array() };
    apply_level_shift(input_array);
    apply_DCT(input_array);

    Q_Table q_table{ gen_wiki_q_table() };

    DU_Array<double> expected_result{{
        {-26.,  -3.,  -6.,   2.,   2.,  -1.,  -0.,   0.},
        {  0.,  -2.,  -4.,   1.,   1.,  -0.,  -0.,   0.},
        { -3.,   1.,   5.,  -1.,  -1.,   0.,   0.,  -0.},
        { -3.,   1.,   2.,  -1.,  -0.,   0.,   0.,   0.},
        {  1.,  -0.,  -0.,  -0.,  -0.,   0.,  -0.,   0.},
        { -0.,   0.,   0.,  -0.,  -0.,   0.,   0.,   0.},
        { -0.,   0.,   0.,  -0.,  -0.,  -0.,   0.,  -0.},
        { -0.,   0.,  -0.,  -0.,  -0.,  -0.,   0.,   0.}
    }};

    // Now perform the quantization
    apply_quantization(input_array, q_table);

    // Verify quantization was performed correctly
    for (size_t ind = 0; ind < input_array.size(); ind++)
    {
        REQUIRE_THAT( input_array[ind], WithinRel(expected_result[ind]) );
    }
}

TEST_CASE( "apply_quantization()::correctly quantizes multiple data units", "[apply_quantization()]" ) {
    DU_Array<double> input_array{ gen_2_DU_Array() };
    apply_level_shift(input_array);
    apply_DCT(input_array);

    Q_Table q_table{ gen_wiki_q_table() };

    DU_Array<double> expected_result{
        {
            {-26.,  -3.,  -6.,   2.,   2.,  -1.,  -0.,   0.},
            {  0.,  -2.,  -4.,   1.,   1.,  -0.,  -0.,   0.},
            { -3.,   1.,   5.,  -1.,  -1.,   0.,   0.,  -0.},
            { -3.,   1.,   2.,  -1.,  -0.,   0.,   0.,   0.},
            {  1.,  -0.,  -0.,  -0.,  -0.,   0.,  -0.,   0.},
            { -0.,   0.,   0.,  -0.,  -0.,   0.,   0.,   0.},
            { -0.,   0.,   0.,  -0.,  -0.,  -0.,   0.,  -0.},
            { -0.,   0.,  -0.,  -0.,  -0.,  -0.,   0.,   0.}
        },
        {
            {-26.,   0.,  -5.,  -3.,   1.,  -0.,  -0.,  -0.},
            { -3.,  -2.,   1.,   1.,  -0.,   0.,   0.,   0.},
            { -4.,  -5.,   5.,   1.,  -0.,   0.,   0.,  -0.},
            {  2.,   1.,  -1.,  -1.,  -0.,  -0.,  -0.,  -0.},
            {  3.,   1.,  -1.,  -0.,  -0.,  -0.,  -0.,  -0.},
            { -1.,  -0.,   0.,   0.,   0.,   0.,  -0.,  -0.},
            { -0.,  -0.,   0.,   0.,  -0.,   0.,   0.,   0.},
            {  0.,   0.,  -0.,   0.,   0.,   0.,  -0.,   0.}
        }
    };

    // Now perform the quantization
    apply_quantization(input_array, q_table);

    // Verify quantization was performed correctly
    for (size_t ind = 0; ind < input_array.size(); ind++)
    {
        REQUIRE_THAT( input_array[ind], WithinRel(expected_result[ind]) );
    }
}

TEST_CASE( "compute_ssss()::computes the ssss value of a number", "[compute_ssss()]" ) {

    SECTION( "n = 0")
    {
        unsigned int n{ 0 };
        unsigned int expected_result{ 0 };
        unsigned int actual_result{ compute_ssss(n) };

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "n = 1")
    {
        unsigned int n{ 1 };
        unsigned int expected_result{ 1 };
        unsigned int actual_result{ compute_ssss(n) };

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "n = 4")
    {
        unsigned int n{ 4 };
        unsigned int expected_result{ 3 };
        unsigned int actual_result{ compute_ssss(n) };

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "n = 5")
    {
        unsigned int n{ 5 };
        unsigned int expected_result{ 3 };
        unsigned int actual_result{ compute_ssss(n) };

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "n = 7")
    {
        unsigned int n{ 7 };
        unsigned int expected_result{ 3 };
        unsigned int actual_result{ compute_ssss(n) };

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "n = 1024")
    {
        unsigned int n{ 1024 };
        unsigned int expected_result{ 11 };
        unsigned int actual_result{ compute_ssss(n) };

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "n = 2047")
    {
        unsigned int n{ 2047 };
        unsigned int expected_result{ 11 };
        unsigned int actual_result{ compute_ssss(n) };

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "n = 2048")
    {
        // Expect this to throw an invalid argument exception
        unsigned int n{ 2048 };

        REQUIRE_THROWS_AS( compute_ssss(n), std::invalid_argument );
    }
}

TEST_CASE( "encode_DC_coeff()::encodes a DC coefficient correctly", "[encode_DC_coeff()]" ) {
    Huff_Table huff_table{ Huff_Table::load_DC_table(Image_Component::Luminance) };
    Bit_String actual_result{};

    SECTION( "diff = 0")
    {
        int diff{ 0 };
        Bit_String expected_result{ "00" };
        
        encode_DC_coeff(actual_result, diff, huff_table);

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "diff = 1")
    {
        int diff{ 1 };
        Bit_String expected_result{ "0101" };
        
        encode_DC_coeff(actual_result, diff, huff_table);

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "diff = -1")
    {
        int diff{ -1 };
        Bit_String expected_result{ "0100" };
        
        encode_DC_coeff(actual_result, diff, huff_table);

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "diff = 2")
    {
        int diff{ 2 };
        Bit_String expected_result{ "01110" };
        
        encode_DC_coeff(actual_result, diff, huff_table);

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "diff = -2")
    {
        int diff{ -2 };
        Bit_String expected_result{ "01101" };
        
        encode_DC_coeff(actual_result, diff, huff_table);

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "diff = 5")
    {
        int diff{ 5 };
        Bit_String expected_result{ "100101" };
        
        encode_DC_coeff(actual_result, diff, huff_table);

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "diff = -5")
    {
        int diff{ -5 };
        Bit_String expected_result{ "100010" };
        
        encode_DC_coeff(actual_result, diff, huff_table);

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "diff = 2047")
    {
        int diff{ 2047 };
        Bit_String expected_result{ "11111111011111111111" };
        
        encode_DC_coeff(actual_result, diff, huff_table);

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "diff = -2047")
    {
        int diff{ -2047 };
        Bit_String expected_result{ "11111111000000000000" };
        
        encode_DC_coeff(actual_result, diff, huff_table);

        REQUIRE( actual_result == expected_result );
    }
}

TEST_CASE( "encode_AC_coeff()::correctly encodes for rrrr = 0", "[encode_AC_coeff()]" ) {
    Huff_Table huff_table{ Huff_Table::load_AC_table(Image_Component::Luminance) };
    Bit_String actual_result{};
    unsigned int rrrr{ 0 };

    SECTION( "rrrr = 0, coeff = 1")
    {
        int coeff{ 1 };
        Bit_String expected_result{ "001" };
        
        encode_AC_coeff(actual_result, coeff, rrrr, huff_table);

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "rrrr = 0, coeff = -1")
    {
        int coeff{ -1 };
        Bit_String expected_result{ "000" };
        
        encode_AC_coeff(actual_result, coeff, rrrr, huff_table);

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "rrrr = 0, coeff = 2")
    {
        int coeff{ 2 };
        Bit_String expected_result{ "0110" };
        
        encode_AC_coeff(actual_result, coeff, rrrr, huff_table);

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "rrrr = 0, coeff = -2")
    {
        int coeff{ -2 };
        Bit_String expected_result{ "0101" };
        
        encode_AC_coeff(actual_result, coeff, rrrr, huff_table);

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "rrrr = 0, coeff = 5")
    {
        int coeff{ 5 };
        Bit_String expected_result{ "100101" };
        
        encode_AC_coeff(actual_result, coeff, rrrr, huff_table);

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "rrrr = 0, coeff = -5")
    {
        int coeff{ -5 };
        Bit_String expected_result{ "100010" };
        
        encode_AC_coeff(actual_result, coeff, rrrr, huff_table);

        REQUIRE( actual_result == expected_result );
    }

    // Note maximum ssss for AC coefficients is 10, unlike for DC coefficients where it's 11
    SECTION( "rrrr = 0, coeff = 1023")
    {
        int coeff{ 1023 };
        Bit_String expected_result{ "11111111100000111111111111" };
        
        encode_AC_coeff(actual_result, coeff, rrrr, huff_table);

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "rrrr = 0, coeff = -1023")
    {
        int coeff{ -1023 };
        Bit_String expected_result{ "11111111100000110000000000" };
        
        encode_AC_coeff(actual_result, coeff, rrrr, huff_table);

        REQUIRE( actual_result == expected_result );
    }
}

TEST_CASE( "encode_AC_coeff()::correctly encodes for rrrr = 5", "[encode_AC_coeff()]" ) {
    Huff_Table huff_table{ Huff_Table::load_AC_table(Image_Component::Luminance) };
    Bit_String actual_result{};
    unsigned int rrrr{ 5 };

    SECTION( "rrrr = 5, coeff = 1")
    {
        int coeff{ 1 };
        Bit_String expected_result{ "11110101" };
        
        encode_AC_coeff(actual_result, coeff, rrrr, huff_table);

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "rrrr = 5, coeff = -1")
    {
        int coeff{ -1 };
        Bit_String expected_result{ "11110100" };
        
        encode_AC_coeff(actual_result, coeff, rrrr, huff_table);

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "rrrr = 5, coeff = 2")
    {
        int coeff{ 2 };
        Bit_String expected_result{ "1111111011110" };
        
        encode_AC_coeff(actual_result, coeff, rrrr, huff_table);

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "rrrr = 5, coeff = -2")
    {
        int coeff{ -2 };
        Bit_String expected_result{ "1111111011101" };
        
        encode_AC_coeff(actual_result, coeff, rrrr, huff_table);

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "rrrr = 5, coeff = 5")
    {
        int coeff{ 5 };
        Bit_String expected_result{ "1111111110011110101" };
        
        encode_AC_coeff(actual_result, coeff, rrrr, huff_table);

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "rrrr = 5, coeff = -5")
    {
        int coeff{ -5 };
        Bit_String expected_result{ "1111111110011110010" };
        
        encode_AC_coeff(actual_result, coeff, rrrr, huff_table);

        REQUIRE( actual_result == expected_result );
    }

    // Note maximum ssss for AC coefficients is 10, unlike for DC coefficients where it's 11
    SECTION( "rrrr = 5, coeff = 1023")
    {
        int coeff{ 1023 };
        Bit_String expected_result{ "11111111101001011111111111" };
        
        encode_AC_coeff(actual_result, coeff, rrrr, huff_table);

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "rrrr = 5, coeff = -1023")
    {
        int coeff{ -1023 };
        Bit_String expected_result{ "11111111101001010000000000" };
        
        encode_AC_coeff(actual_result, coeff, rrrr, huff_table);

        REQUIRE( actual_result == expected_result );
    }
}

TEST_CASE( "encode_AC_coeff()::correctly encodes for rrrr = 15", "[encode_AC_coeff()]" ) {
    Huff_Table huff_table{ Huff_Table::load_AC_table(Image_Component::Luminance) };
    Bit_String actual_result{};
    unsigned int rrrr{ 15 };

    SECTION( "rrrr = 15, coeff = 1")
    {
        int coeff{ 1 };
        Bit_String expected_result{ "11111111111101011" };
        
        encode_AC_coeff(actual_result, coeff, rrrr, huff_table);

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "rrrr = 15, coeff = -1")
    {
        int coeff{ -1 };
        Bit_String expected_result{ "11111111111101010" };
        
        encode_AC_coeff(actual_result, coeff, rrrr, huff_table);

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "rrrr = 15, coeff = 2")
    {
        int coeff{ 2 };
        Bit_String expected_result{ "111111111111011010" };
        
        encode_AC_coeff(actual_result, coeff, rrrr, huff_table);

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "rrrr = 15, coeff = -2")
    {
        int coeff{ -2 };
        Bit_String expected_result{ "111111111111011001" };
        
        encode_AC_coeff(actual_result, coeff, rrrr, huff_table);

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "rrrr = 15, coeff = 5")
    {
        int coeff{ 5 };
        Bit_String expected_result{ "1111111111110111101" };
        
        encode_AC_coeff(actual_result, coeff, rrrr, huff_table);

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "rrrr = 15, coeff = -5")
    {
        int coeff{ -5 };
        Bit_String expected_result{ "1111111111110111010" };
        
        encode_AC_coeff(actual_result, coeff, rrrr, huff_table);

        REQUIRE( actual_result == expected_result );
    }

    // Note maximum ssss for AC coefficients is 10, unlike for DC coefficients where it's 11
    SECTION( "rrrr = 15, coeff = 1023")
    {
        int coeff{ 1023 };
        Bit_String expected_result{ "11111111111111101111111111" };
        
        encode_AC_coeff(actual_result, coeff, rrrr, huff_table);

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "rrrr = 15, coeff = -1023")
    {
        int coeff{ -1023 };
        Bit_String expected_result{ "11111111111111100000000000" };
        
        encode_AC_coeff(actual_result, coeff, rrrr, huff_table);

        REQUIRE( actual_result == expected_result );
    }
}

TEST_CASE( "encode_AC_coeffs()::all entries the same", "[encode_AC_coeffs()]" ) {
    Huff_Table huff_table{ Huff_Table::load_AC_table(Image_Component::Luminance) };
    size_t du_ind{ 0 };
    Bit_String actual_result{};

    SECTION( "all entries zero") {
        DU_Array<double> input_array{{
            {0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 0, 0, 0},
        }};
        
        Bit_String expected_result{"1010"};

        encode_AC_coeffs(actual_result, input_array, du_ind, huff_table);

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "all entries 1") {
        DU_Array<double> input_array{{
            {1, 1, 1, 1, 1, 1, 1, 1},
            {1, 1, 1, 1, 1, 1, 1, 1},
            {1, 1, 1, 1, 1, 1, 1, 1},
            {1, 1, 1, 1, 1, 1, 1, 1},
            {1, 1, 1, 1, 1, 1, 1, 1},
            {1, 1, 1, 1, 1, 1, 1, 1},
            {1, 1, 1, 1, 1, 1, 1, 1},
            {1, 1, 1, 1, 1, 1, 1, 1},
        }};
        // Expected to be 001 repeated 63 times
        Bit_String expected_result{
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
            "001"
        };

        encode_AC_coeffs(actual_result, input_array, du_ind, huff_table);

        REQUIRE( actual_result == expected_result );
    }
    SECTION( "all entries 1") {
        DU_Array<double> input_array{{
            {-1, -1, -1, -1, -1, -1, -1, -1},
            {-1, -1, -1, -1, -1, -1, -1, -1},
            {-1, -1, -1, -1, -1, -1, -1, -1},
            {-1, -1, -1, -1, -1, -1, -1, -1},
            {-1, -1, -1, -1, -1, -1, -1, -1},
            {-1, -1, -1, -1, -1, -1, -1, -1},
            {-1, -1, -1, -1, -1, -1, -1, -1},
            {-1, -1, -1, -1, -1, -1, -1, -1},
        }};
        // Expected to be 000 repeated 63 times
        Bit_String expected_result{
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
            "000"
        };

        encode_AC_coeffs(actual_result, input_array, du_ind, huff_table);

        REQUIRE( actual_result == expected_result );
    }
}

TEST_CASE( "encode_AC_coeffs()::zig-zag ordering", "[encode_AC_coeffs()]" ) {
    Huff_Table huff_table{ Huff_Table::load_AC_table(Image_Component::Luminance) };
    size_t du_ind{ 0 };
    
    DU_Array<double> input_array{{
        {0, 1, 5, 0, 0, 0, 0, 0},
        {2, 4, 0, 0, 0, 0, 0, 0},
        {3, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
    }};

    Bit_String expected_result{
        "001"     // 1
        "0110"    // 2
        "0111"    // 3
        "100100"  // 4
        "100101"  // 5
        "1010"    // EOB
    };
    
    Bit_String actual_result{};

    encode_AC_coeffs(actual_result, input_array, du_ind, huff_table);

    REQUIRE( actual_result == expected_result );
}

TEST_CASE( "encode_AC_coeffs()::run of 1 zero", "[encode_AC_coeffs()]" ) {
    Huff_Table huff_table{ Huff_Table::load_AC_table(Image_Component::Luminance) };
    size_t du_ind{ 0 };
    
    DU_Array<double> input_array{{
        {0, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {1, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
    }};

    Bit_String expected_result{
        "001"    // 1
        "11001"  // 1, with run of 1 zero
        "1010"   // EOB
    };
    
    Bit_String actual_result{};

    encode_AC_coeffs(actual_result, input_array, du_ind, huff_table);

    REQUIRE( actual_result == expected_result );
}

TEST_CASE( "encode_AC_coeffs()::run of 5 zeros", "[encode_AC_coeffs()]" ) {
    Huff_Table huff_table{ Huff_Table::load_AC_table(Image_Component::Luminance) };
    size_t du_ind{ 0 };
    
    DU_Array<double> input_array{{
        {0, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
    }};

    Bit_String expected_result{
        "001"       // 1
        "11110101"  // 1, with run of 5 zeros
        "1010"      // EOB
    };
    
    Bit_String actual_result{};

    encode_AC_coeffs(actual_result, input_array, du_ind, huff_table);

    REQUIRE( actual_result == expected_result );
}

TEST_CASE( "encode_AC_coeffs()::run of 15 zeros", "[encode_AC_coeffs()]" ) {
    Huff_Table huff_table{ Huff_Table::load_AC_table(Image_Component::Luminance) };
    size_t du_ind{ 0 };
    
    DU_Array<double> input_array{{
        {0, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
    }};

    Bit_String expected_result{
        "001"                // 1
        "11111111111101011"  // 1, with run of 15 zeros
        "1010"               // EOB
    };
    
    Bit_String actual_result{};

    encode_AC_coeffs(actual_result, input_array, du_ind, huff_table);

    REQUIRE( actual_result == expected_result );
}

TEST_CASE( "encode_AC_coeffs()::run of 16 zeros", "[encode_AC_coeffs()]" ) {
    Huff_Table huff_table{ Huff_Table::load_AC_table(Image_Component::Luminance) };
    size_t du_ind{ 0 };
    
    DU_Array<double> input_array{{
        {0, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
    }};

    Bit_String expected_result{
        "001"          // 1
        "11111111001"  // ZRL (run of 16 zeros)
        "001"          // 1
        "1010"         // EOB
    };
    
    Bit_String actual_result{};

    encode_AC_coeffs(actual_result, input_array, du_ind, huff_table);

    REQUIRE( actual_result == expected_result );
}

TEST_CASE( "encode_AC_coeffs()::run of 17 zeros", "[encode_AC_coeffs()]" ) {
    Huff_Table huff_table{ Huff_Table::load_AC_table(Image_Component::Luminance) };
    size_t du_ind{ 0 };
    
    DU_Array<double> input_array{{
        {0, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
    }};

    Bit_String expected_result{
        "001"          // 1
        "11111111001"  // ZRL (run of 16 zeros)
        "11001"        // 1, with run of 1 zero
        "1010"         // EOB
    };
    
    Bit_String actual_result{};

    encode_AC_coeffs(actual_result, input_array, du_ind, huff_table);

    REQUIRE( actual_result == expected_result );
}

TEST_CASE( "encode_AC_coeffs()::run of 32 zeros", "[encode_AC_coeffs()]" ) {
    Huff_Table huff_table{ Huff_Table::load_AC_table(Image_Component::Luminance) };
    size_t du_ind{ 0 };
    
    DU_Array<double> input_array{{
        {0, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0},
    }};

    Bit_String expected_result{
        "001"          // 1
        "11111111001"  // ZRL (run of 16 zeros)
        "11111111001"  // ZRL (run of 16 zeros)
        "001"          // 1
        "1010"         // EOB
    };
    
    Bit_String actual_result{};

    encode_AC_coeffs(actual_result, input_array, du_ind, huff_table);

    REQUIRE( actual_result == expected_result );
}

TEST_CASE( "encode_data_unit_sequential()::encode sample data unit", "[encode_data_unit_sequential()]" ) {
    Huff_Table huff_dc{ Huff_Table::load_DC_table(Image_Component::Luminance) };
    Huff_Table huff_ac{ Huff_Table::load_AC_table(Image_Component::Luminance) };
    size_t du_ind{ 0 };
    int prev_dc{ 0 };

    // Corresponds to the Wikipedia example using the Wikipedia quantization table
    DU_Array<double> input_array{{
        {-26.,  -3.,  -6.,   2.,   2.,  -1.,  -0.,   0.},
        {  0.,  -2.,  -4.,   1.,   1.,  -0.,  -0.,   0.},
        { -3.,   1.,   5.,  -1.,  -1.,   0.,   0.,  -0.},
        { -3.,   1.,   2.,  -1.,  -0.,   0.,   0.,   0.},
        {  1.,  -0.,  -0.,  -0.,  -0.,   0.,  -0.,   0.},
        { -0.,   0.,   0.,  -0.,  -0.,   0.,   0.,   0.},
        { -0.,   0.,   0.,  -0.,  -0.,  -0.,   0.,  -0.},
        { -0.,   0.,  -0.,  -0.,  -0.,  -0.,   0.,   0.}
    }};

    Bit_String expected_result{
        "11000101"  // -26, ssss=5
        "0100"      // -3, rs = 0, 2
        "1101100"   // -3, rs = 1, 2
        "0101"      // -2  rs = 0, 2
        "100001"    // -6  rs = 0, 3
        "0110"      //  2  rs = 0, 2
        "100011"    // -4  rs = 0, 3
        "001"       //  1  rs = 0, 1
        "0100"      // -3  rs = 0, 2
        "001"       //  1  rs = 0, 1
        "001"       //  1  rs = 0, 1
        "100101"    //  5  rs = 0, 3
        "001"       //  1  rs = 0, 1
        "0110"      //  2  rs = 0, 2
        "000"       // -1  rs = 0, 1
        "001"       //  1  rs = 0, 1
        "000"       // -1  rs = 0, 1
        "0110"      //  2  rs = 0, 2
        "11110100"  // -1  rs = 5, 1
        "000"       // -1  rs = 0, 1
        "1010"      // EOB
    };

    Bit_String actual_result{};

    encode_data_unit_sequential(actual_result, input_array, du_ind, prev_dc, huff_dc, huff_ac);

    REQUIRE( actual_result == expected_result );
}

TEST_CASE( "append_q_table_marker_segment()::single quantization table", "[append_q_table_marker_segment()]" ) {
    std::vector<Q_Table> q_tables{ gen_wiki_q_table() };
    std::vector<unsigned int> destination_indices{ 0 };

    std::vector<unsigned char> expected_result = {
        0xFF,   // DQT segment marker
        0xDB,
        0x00,   // Length of marker segment (2 bytes for length, 1 for table type/destination
        67,     // and 64 for the table's elements)
        0x00,   //  Pq/Tq byte, 8 bit precision & destination identifier
        16,     // Quantization table elements in zig-zag order
        11,
        12,
        14,
        12,
        10,
        16,
        14,
        13,
        14,
        18,
        17,
        16,
        19,
        24,
        40,
        26,
        24,
        22,
        22,
        24,
        49,
        35,
        37,
        29,
        40,
        58,
        51,
        61,
        60,
        57,
        51,
        56,
        55,
        64,
        72,
        92,
        78,
        64,
        68,
        87,
        69,
        55,
        56,
        80,
        109,
        81,
        87,
        95,
        98,
        103,
        104,
        103,
        62,
        77,
        113,
        121,
        112,
        100,
        120,
        92,
        101,
        103,
        99
    };

    std::vector<unsigned char> actual_result{};
    append_q_table_marker_segment(actual_result, q_tables, destination_indices);

    REQUIRE_THAT( actual_result, RangeEquals(expected_result) );
}

TEST_CASE( "append_q_table_marker_segment()::two quantization table", "[append_q_table_marker_segment()]" ) {
    Q_Table wiki_q_table{ gen_wiki_q_table() };
    std::vector<Q_Table> q_tables(2);

    // Create a vector of q tables where the first is the wiki table, the second is
    // is derived from the wiki table with +1 added elementwise
    for (size_t table_ind = 0; table_ind < q_tables.size(); table_ind++)
    {
        for (size_t i = 0; i < q_tables[table_ind].shape()[1]; i++)
        {
            for (size_t j = 0; j < q_tables[table_ind].shape()[1]; j++)
            {
                q_tables[table_ind](i, j) = wiki_q_table(i, j) + table_ind;
            }
        }
    }
    
    std::vector<unsigned int> destination_indices{ 0, 1 };

    std::vector<unsigned char> expected_result = {
        0xFF,   // DQT segment marker
        0xDB,
        0x00,   // Length of marker segment (2 bytes for length, 2 for table type/destination
        132,     // and 128 for the table's elements)
        0x00,   // Pq/Tq byte, 8 bit precision & destination identifier
        16,     // Quantization table elements in zig-zag order
        11,
        12,
        14,
        12,
        10,
        16,
        14,
        13,
        14,
        18,
        17,
        16,
        19,
        24,
        40,
        26,
        24,
        22,
        22,
        24,
        49,
        35,
        37,
        29,
        40,
        58,
        51,
        61,
        60,
        57,
        51,
        56,
        55,
        64,
        72,
        92,
        78,
        64,
        68,
        87,
        69,
        55,
        56,
        80,
        109,
        81,
        87,
        95,
        98,
        103,
        104,
        103,
        62,
        77,
        113,
        121,
        112,
        100,
        120,
        92,
        101,
        103,
        99,     // Last element of first q table
        0x01,   // Pq/Tq byte, 8 bit precision & destination identifier
        17,     // Quantization table elements in zig-zag order
        12,
        13,
        15,
        13,
        11,
        17,
        15,
        14,
        15,
        19,
        18,
        17,
        20,
        25,
        41,
        27,
        25,
        23,
        23,
        25,
        50,
        36,
        38,
        30,
        41,
        59,
        52,
        62,
        61,
        58,
        52,
        57,
        56,
        65,
        73,
        93,
        79,
        65,
        69,
        88,
        70,
        56,
        57,
        81,
        110,
        82,
        88,
        96,
        99,
        104,
        105,
        104,
        63,
        78,
        114,
        122,
        113,
        101,
        121,
        93,
        102,
        104,
        100
    };

    std::vector<unsigned char> actual_result{};
    append_q_table_marker_segment(actual_result, q_tables, destination_indices);

    REQUIRE_THAT( actual_result, RangeEquals(expected_result) );
}

TEST_CASE( "append_huff_table_marker_segment()::one Huffman table", "[append_huff_table_marker_segment()]" ) {
    Huff_Table luminance_table{ Huff_Table::load_DC_table(Image_Component::Luminance) };
    std::vector<Huff_Table_Ref> huff_tables{
        {luminance_table, Huff_Table_Ref::Huff_Table_Type::DC, 0}
    };

    std::vector<unsigned char> expected_result{
        0xFF,   // DHT marker, FFC4
        0xC4,
        0,      // Length, 2 bytes, most significant first
        31,
        0x00,   // Table class (0) and destination (0)
        0,      // Start of BITS array
        1,
        5,
        1,
        1,
        1,
        1,
        1,
        1,
        0,
        0,
        0,
        0,
        0,
        0,
        0,      // End of BITS array
        0,      // Start of HUFFVAL array
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11      // End of HUFFVAL array 
    };

    std::vector<unsigned char> actual_result;

    append_huff_table_marker_segment(actual_result, huff_tables);

    REQUIRE_THAT( actual_result, RangeEquals(expected_result) );
}

TEST_CASE( "append_huff_table_marker_segment()::two Huffman tables", "[append_huff_table_marker_segment()]" ) {
    Huff_Table luminance_table{ Huff_Table::load_DC_table(Image_Component::Luminance) };
    Huff_Table chromiance_table{ Huff_Table::load_DC_table(Image_Component::Chrominance) };
    std::vector<Huff_Table_Ref> huff_tables{
        {luminance_table, Huff_Table_Ref::Huff_Table_Type::DC, 0},
        // We'll pretend the DC chromiance table is AC to make the test shorter
        {chromiance_table, Huff_Table_Ref::Huff_Table_Type::AC, 1},
    };

    std::vector<unsigned char> expected_result{
        0xFF,   // DHT marker, FFC4
        0xC4,
        0,      // Length, 2 bytes, most significant first
        60,
        // Start of first table data
        0x00,   // Table class (0) and destination (0)
        0,      // Start of BITS array
        1,
        5,
        1,
        1,
        1,
        1,
        1,
        1,
        0,
        0,
        0,
        0,
        0,
        0,
        0,      // End of BITS array
        0,      // Start of HUFFVAL array
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,     // End of HUFFVAL array
        // Start of second table data
        0x11,   // Table class (0) and destination (0)
        0,      // Start of BITS array
        3,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        0,
        0,
        0,
        0,
        0,      // End of BITS array
        0,      // Start of HUFFVAL array
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,     // End of HUFFVAL array
    };

    std::vector<unsigned char> actual_result;

    append_huff_table_marker_segment(actual_result, huff_tables);

    REQUIRE_THAT( actual_result, RangeEquals(expected_result) );
}

TEST_CASE( "append_mcu()::single component", "[append_mcu()]" ) {
    std::vector<int> prev_dc{ 0 };
    std::vector<size_t> du_ind{ 0 };
    std::vector<Comp_Info> comp_info{
        Comp_Info(
            0,  // Q table ind
            0,  // DC Huffman ind
            0,  // AC Huffman ind
            1,  // H sampling factor
            1   // V sampling factor
        )
    };
    std::vector<Huff_Table> huff_dc{ Huff_Table::load_DC_table(Image_Component::Luminance) };
    std::vector<Huff_Table> huff_ac{ Huff_Table::load_AC_table(Image_Component::Luminance) };

    // Corresponds to the Wikipedia example using the Wikipedia quantization table
    std::vector<DU_Array<double>> input_array{{
        {
            {-26.,  -3.,  -6.,   2.,   2.,  -1.,  -0.,   0.},
            {  0.,  -2.,  -4.,   1.,   1.,  -0.,  -0.,   0.},
            { -3.,   1.,   5.,  -1.,  -1.,   0.,   0.,  -0.},
            { -3.,   1.,   2.,  -1.,  -0.,   0.,   0.,   0.},
            {  1.,  -0.,  -0.,  -0.,  -0.,   0.,  -0.,   0.},
            { -0.,   0.,   0.,  -0.,  -0.,   0.,   0.,   0.},
            { -0.,   0.,   0.,  -0.,  -0.,  -0.,   0.,  -0.},
            { -0.,   0.,  -0.,  -0.,  -0.,  -0.,   0.,   0.}
        },
        {
            {-26.,  -3.,  -6.,   2.,   2.,  -1.,  -0.,   0.},
            {  0.,  -2.,  -4.,   1.,   1.,  -0.,  -0.,   0.},
            { -3.,   1.,   5.,  -1.,  -1.,   0.,   0.,  -0.},
            { -3.,   1.,   2.,  -1.,  -0.,   0.,   0.,   0.},
            {  1.,  -0.,  -0.,  -0.,  -0.,   0.,  -0.,   0.},
            { -0.,   0.,   0.,  -0.,  -0.,   0.,   0.,   0.},
            { -0.,   0.,   0.,  -0.,  -0.,  -0.,   0.,  -0.},
            { -0.,   0.,  -0.,  -0.,  -0.,  -0.,   0.,   0.}
        }

    }};

    Bit_String expected_result{
        // First data unit
        "11000101"  // -26, ssss=5
        "0100"      // -3, rs = 0, 2
        "1101100"   // -3, rs = 1, 2
        "0101"      // -2  rs = 0, 2
        "100001"    // -6  rs = 0, 3
        "0110"      //  2  rs = 0, 2
        "100011"    // -4  rs = 0, 3
        "001"       //  1  rs = 0, 1
        "0100"      // -3  rs = 0, 2
        "001"       //  1  rs = 0, 1
        "001"       //  1  rs = 0, 1
        "100101"    //  5  rs = 0, 3
        "001"       //  1  rs = 0, 1
        "0110"      //  2  rs = 0, 2
        "000"       // -1  rs = 0, 1
        "001"       //  1  rs = 0, 1
        "000"       // -1  rs = 0, 1
        "0110"      //  2  rs = 0, 2
        "11110100"  // -1  rs = 5, 1
        "000"       // -1  rs = 0, 1
        "1010"      // EOB
        // Second data unit
        // Remember the DC elements of both data units are the same and its
        // the difference that is encoded. I.e. 0
        "00"        // 0, ssss=5
        "0100"      // -3, rs = 0, 2
        "1101100"   // -3, rs = 1, 2
        "0101"      // -2  rs = 0, 2
        "100001"    // -6  rs = 0, 3
        "0110"      //  2  rs = 0, 2
        "100011"    // -4  rs = 0, 3
        "001"       //  1  rs = 0, 1
        "0100"      // -3  rs = 0, 2
        "001"       //  1  rs = 0, 1
        "001"       //  1  rs = 0, 1
        "100101"    //  5  rs = 0, 3
        "001"       //  1  rs = 0, 1
        "0110"      //  2  rs = 0, 2
        "000"       // -1  rs = 0, 1
        "001"       //  1  rs = 0, 1
        "000"       // -1  rs = 0, 1
        "0110"      //  2  rs = 0, 2
        "11110100"  // -1  rs = 5, 1
        "000"       // -1  rs = 0, 1
        "1010"      // EOB
    };

    Bit_String actual_result{};
    size_t mcu_count{ 2 };

    for (size_t mcu_ind = 0; mcu_ind < mcu_count; mcu_ind++)
    {
        append_mcu(actual_result, prev_dc, du_ind, input_array, comp_info, huff_dc, huff_ac);
    }

    REQUIRE( actual_result == expected_result );
}

TEST_CASE( "append_mcu()::multiple components", "[append_mcu()]" ) {
    std::vector<int> prev_dc{ 0, 0, 0 };
    std::vector<size_t> du_ind{ 0, 0, 0 };
    std::vector<Comp_Info> comp_info{
        Comp_Info(
            0,  // Q table ind
            0,  // DC Huffman ind
            0,  // AC Huffman ind
            2,  // H sampling factor
            2   // V sampling factor
        ),
        Comp_Info(
            0,  // Q table ind
            1,  // DC Huffman ind
            1,  // AC Huffman ind
            1,  // H sampling factor
            1   // V sampling factor
        ),
        Comp_Info(
            0,  // Q table ind
            1,  // DC Huffman ind
            1,  // AC Huffman ind
            1,  // H sampling factor
            1   // V sampling factor
        )
    };
    std::vector<Huff_Table> huff_dc{ 
        Huff_Table::load_DC_table(Image_Component::Luminance),
        Huff_Table::load_DC_table(Image_Component::Chrominance)
    };
    std::vector<Huff_Table> huff_ac{ 
        Huff_Table::load_AC_table(Image_Component::Luminance),
        Huff_Table::load_AC_table(Image_Component::Chrominance)
    };

    // Corresponds to the Wikipedia example using the Wikipedia quantization table
    std::vector<DU_Array<double>> input_array{
        {
            // First component has 4 data units
            {
                { 1, 2, 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0, 0, 0 },
            },
            {
                { 3, 4, 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0, 0, 0 },
            },
            {
                { 5, 6, 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0, 0, 0 },
            },
            {
                { 7, 8, 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0, 0, 0 },
                { 0, 0, 0, 0, 0, 0, 0, 0 },
            }
        },
        // Next 2 components only have 1 data unit each
        {
            {
                { 9, 10, 0, 0, 0, 0, 0, 0 },
                { 0,  0, 0, 0, 0, 0, 0, 0 },
                { 0,  0, 0, 0, 0, 0, 0, 0 },
                { 0,  0, 0, 0, 0, 0, 0, 0 },
                { 0,  0, 0, 0, 0, 0, 0, 0 },
                { 0,  0, 0, 0, 0, 0, 0, 0 },
                { 0,  0, 0, 0, 0, 0, 0, 0 },
                { 0,  0, 0, 0, 0, 0, 0, 0 },
            }
        },
        {
            {
                { 11, 12, 0, 0, 0, 0, 0, 0 },
                {  0,  0, 0, 0, 0, 0, 0, 0 },
                {  0,  0, 0, 0, 0, 0, 0, 0 },
                {  0,  0, 0, 0, 0, 0, 0, 0 },
                {  0,  0, 0, 0, 0, 0, 0, 0 },
                {  0,  0, 0, 0, 0, 0, 0, 0 },
                {  0,  0, 0, 0, 0, 0, 0, 0 },
                {  0,  0, 0, 0, 0, 0, 0, 0 },
            }
        }
    };

    Bit_String expected_result{
        // First component
        // First data unit
        "0101"      // 1, ssss=1
        "0110"      // 2, rs = 0, 2
        "1010"      // EOB
        // Second data unit
        "01110"     // 2, ssss=2  (note diff=2)
        "100100"    // 4, rs = 0, 3
        "1010"      // EOB
        // Third data unit
        "01110"     // 2, ssss=3
        "100110"    // 6, rs = 0, 3
        "1010"      // EOB
        // Fourth data unit
        "01110"     // 2, ssss=3
        "10111000"  // 8, rs = 0, 4
        "1010"      // EOB

        // Second component
        "11101001"  // 9, ssss=4
        "110001010" // 10, rs = 0, 4
        "00"      // EOB

        // Third component
        "11101011"  // 11, ssss=4
        "110001100" // 12, rs = 0, 4
        "00"      // EOB
    };

    Bit_String actual_result{};

    // Append the MCU
    append_mcu(actual_result, prev_dc, du_ind, input_array, comp_info, huff_dc, huff_ac);

    REQUIRE( actual_result == expected_result );
}

TEST_CASE( "append_scan_header()::single component", "[append_scan_header()]" ) {

    std::vector<Comp_Info> comp_infos{
        Comp_Info{0, 0, 0, 1, 1}
    };

    std::vector<unsigned char> expected_result{
        0xFF,   // DHT marker, FFC4
        0xDA,
        0,      // Length, 2 bytes, most significant first
        8,
        1,      // Number of components in scan
        // First component
        0,      // Component identifier
        0x00,   // Composite byte of DC table destination (most significant) and AC table destination
        // End of component data
        0,      // Start of spectral or predictor selection (0 for sequential DCT)
        63,     // End of spectral selection (63 for sequential DCT)
        0x00    // Composite byte of successive approximation bit position high/low,
                // 0 in both cases
    };

    std::vector<unsigned char> actual_result;

    append_scan_header(actual_result, comp_infos);

    REQUIRE_THAT( actual_result, RangeEquals(expected_result) );
}

TEST_CASE( "append_scan_header()::multiple components", "[append_scan_header()]" ) {

    std::vector<Comp_Info> comp_infos{
        Comp_Info{0, 0, 0, 2, 1},
        Comp_Info{0, 1, 1, 1, 1},
        Comp_Info{0, 1, 1, 1, 1}
    };

    std::vector<unsigned char> expected_result{
        0xFF,   // DHT marker, FFC4
        0xDA,
        0,      // Length, 2 bytes, most significant first
        12,
        3,      // Number of components in scan
        // First component
        0,      // Component identifier
        0x00,   // Composite byte of DC table destination (most significant) and AC table destination
        // Second component
        1,      // Component identifier
        0x11,   // Composite byte of DC table destination (most significant) and AC table destination
        // Third component
        2,      // Component identifier
        0x11,   // Composite byte of DC table destination (most significant) and AC table destination
        // End of component data
        0,      // Start of spectral or predictor selection (0 for sequential DCT)
        63,     // End of spectral selection (63 for sequential DCT)
        0x00    // Composite byte of successive approximation bit position high/low,
                // 0 in both cases
    };

    std::vector<unsigned char> actual_result;

    append_scan_header(actual_result, comp_infos);

    REQUIRE_THAT( actual_result, RangeEquals(expected_result) );
}
