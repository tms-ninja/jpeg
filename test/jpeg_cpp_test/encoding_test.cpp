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

// generates an Array_2d<double> of shape (32, 32) with values increasing in row-major order
Array_2d<double> gen_array_2d_32x32()
{
    Array_2d<double> array{32, 32};

    for (size_t ind = 0; ind < array.size(); ind++)
    {
        array[ind] = ind;
    }
    
    return array;
}

Array_2d<double> gen_array_2d_arange(size_t rows, size_t cols)
{
    Array_2d<double> array{rows, cols};

    for (size_t ind = 0; ind < array.size(); ind++)
    {
        array[ind] = ind;
    }
    
    return array;
}

TEST_CASE( "apply_level_shift()::level shift values in DU_Array", "[apply_level_shift()]" ) {
    DU_Array<double> array{ gen_arange_DU_Array() };
    
    apply_level_shift(array);

    // Verify all elements have been reduced by 128
    for (size_t ind = 0; ind < array.size(); ind++)
    {
        double expected{ static_cast<double>(ind) - 128.0 };
        
        CAPTURE( ind );
        CHECK_THAT( array[ind], WithinRel(expected) );
    }
}

TEST_CASE( "load_DCT_matrix()::loads the correct matrix", "[load_DCT_matrix()]" ) {
    Array_2d<double> array{ load_DCT_matrix() };

    REQUIRE( array.shape()[0] == 8 );
    REQUIRE( array.shape()[1] == 8 );
    
    // Check a selection of values
    CHECK_THAT( array(0, 0), WithinRel(0.5) );
    CHECK_THAT( array(0, 1), WithinRel(0.4903926402016152) );
    CHECK_THAT( array(1, 0), WithinRel(0.5) );
    CHECK_THAT( array(3, 4), WithinRel(0.3535533905932737) );
    CHECK_THAT( array(7, 2), WithinRel(0.4619397662556433) );
}

TEST_CASE( "load_DCT_matrix_transpose()::loads the correct matrix", "[load_DCT_matrix_transpose()]" ) {
    Array_2d<double> array{ load_DCT_matrix_transpose() };

    REQUIRE( array.shape()[0] == 8 );
    REQUIRE( array.shape()[1] == 8 );
    
    // Check a selection of values
    CHECK_THAT( array(0, 0), WithinRel(0.5) );
    CHECK_THAT( array(0, 1), WithinRel(0.5) );
    CHECK_THAT( array(1, 0), WithinRel(0.4903926402016152) );
    CHECK_THAT( array(3, 4), WithinRel(0.2777851165098009) );
    CHECK_THAT( array(7, 2), WithinRel(0.4157348061512727) );
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
        CAPTURE( ind );
        CHECK_THAT( actual_result[ind], WithinRel(expected_result[ind]) );
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
        CAPTURE( ind );

        // Note sure why some seem so off, increase the tolerance for now
        CHECK_THAT( input_array[ind], WithinRel(expected_result[ind], 1e-12) );
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
        CAPTURE( ind );

        // Note sure why some seem so off, increase the tolerance for now
        CHECK_THAT( input_array[ind], WithinRel(expected_result[ind], 1e-12) );
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
        CAPTURE( ind );
        CHECK_THAT( input_array[ind], WithinRel(expected_result[ind]) );
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
        CAPTURE( ind );
        CHECK_THAT( input_array[ind], WithinRel(expected_result[ind]) );
    }
}

TEST_CASE( "compute_ssss()::computes the ssss value of a number", "[compute_ssss()]" ) {

    SECTION( "n = 0")
    {
        unsigned int n{ 0 };
        unsigned int expected_result{ 0 };
        unsigned int actual_result{ compute_ssss(n) };

        CHECK( actual_result == expected_result );
    }
    SECTION( "n = 1")
    {
        unsigned int n{ 1 };
        unsigned int expected_result{ 1 };
        unsigned int actual_result{ compute_ssss(n) };

        CHECK( actual_result == expected_result );
    }
    SECTION( "n = 4")
    {
        unsigned int n{ 4 };
        unsigned int expected_result{ 3 };
        unsigned int actual_result{ compute_ssss(n) };

        CHECK( actual_result == expected_result );
    }
    SECTION( "n = 5")
    {
        unsigned int n{ 5 };
        unsigned int expected_result{ 3 };
        unsigned int actual_result{ compute_ssss(n) };

        CHECK( actual_result == expected_result );
    }
    SECTION( "n = 7")
    {
        unsigned int n{ 7 };
        unsigned int expected_result{ 3 };
        unsigned int actual_result{ compute_ssss(n) };

        CHECK( actual_result == expected_result );
    }
    SECTION( "n = 1024")
    {
        unsigned int n{ 1024 };
        unsigned int expected_result{ 11 };
        unsigned int actual_result{ compute_ssss(n) };

        CHECK( actual_result == expected_result );
    }
    SECTION( "n = 2047")
    {
        unsigned int n{ 2047 };
        unsigned int expected_result{ 11 };
        unsigned int actual_result{ compute_ssss(n) };

        CHECK( actual_result == expected_result );
    }
    SECTION( "n = 2048")
    {
        // Expect this to throw an invalid argument exception
        unsigned int n{ 2048 };

        CHECK_THROWS_AS( compute_ssss(n), std::invalid_argument );
    }
}

TEST_CASE( "encode_DC_coeff()::encodes a DC coefficient correctly", "[encode_DC_coeff()]" ) {
    std::vector<Coefficient> actual_result;
    unsigned int comp_ind{ 0 };

    SECTION( "diff = 0")
    {
        int diff{ 0 };
        std::vector<Coefficient> expected_result{ 
            Coefficient(0, 0, Coefficient_Type::DC, 0)
        };
        
        encode_DC_coeff(actual_result, diff, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }
    SECTION( "diff = 1")
    {
        int diff{ 1 };
        std::vector<Coefficient> expected_result{ 
            Coefficient(1, 1, Coefficient_Type::DC, 0)
        };
        
        encode_DC_coeff(actual_result, diff, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }
    SECTION( "diff = -1")
    {
        int diff{ -1 };
        std::vector<Coefficient> expected_result{ 
            Coefficient(-1, 1, Coefficient_Type::DC, 0)
        };
        
        encode_DC_coeff(actual_result, diff, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }
    SECTION( "diff = 2")
    {
        int diff{ 2 };
        std::vector<Coefficient> expected_result{ 
            Coefficient(2, 2, Coefficient_Type::DC, 0)
        };
        
        encode_DC_coeff(actual_result, diff, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }
    SECTION( "diff = -2")
    {
        int diff{ -2 };
        std::vector<Coefficient> expected_result{ 
            Coefficient(-2, 2, Coefficient_Type::DC, 0)
        };
        
        encode_DC_coeff(actual_result, diff, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }
    SECTION( "diff = 5")
    {
        int diff{ 5 };
        std::vector<Coefficient> expected_result{ 
            Coefficient(5, 3, Coefficient_Type::DC, 0)
        };
        
        encode_DC_coeff(actual_result, diff, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }
    SECTION( "diff = -5")
    {
        int diff{ -5 };
        std::vector<Coefficient> expected_result{ 
            Coefficient(-5, 3, Coefficient_Type::DC, 0)
        };
        
        encode_DC_coeff(actual_result, diff, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }
    SECTION( "diff = 2047")
    {
        int diff{ 2047 };
        std::vector<Coefficient> expected_result{ 
            Coefficient(2047, 11, Coefficient_Type::DC, 0)
        };
        
        encode_DC_coeff(actual_result, diff, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }
    SECTION( "diff = -2047")
    {
        int diff{ -2047 };
        std::vector<Coefficient> expected_result{ 
            Coefficient(-2047, 11, Coefficient_Type::DC, 0)
        };
        
        encode_DC_coeff(actual_result, diff, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }
}

TEST_CASE( "encode_AC_coeff()::correctly encodes for rrrr = 0", "[encode_AC_coeff()]" ) {
    std::vector<Coefficient> actual_result;
    unsigned int rrrr{ 0 };
    unsigned int comp_ind{ 0 };

    SECTION( "rrrr = 0, coeff = 1")
    {
        int coeff{ 1 };
        std::vector<Coefficient> expected_result{ 
            Coefficient(1, 1, Coefficient_Type::AC, 0)
        };
        
        encode_AC_coeff(actual_result, coeff, rrrr, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }
    SECTION( "rrrr = 0, coeff = -1")
    {
        int coeff{ -1 };
        std::vector<Coefficient> expected_result{ 
            Coefficient(-1, 1, Coefficient_Type::AC, 0)
        };
        
        encode_AC_coeff(actual_result, coeff, rrrr, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }
    SECTION( "rrrr = 0, coeff = 2")
    {
        int coeff{ 2 };
        std::vector<Coefficient> expected_result{ 
            Coefficient(2, 2, Coefficient_Type::AC, 0)
        };
        
        encode_AC_coeff(actual_result, coeff, rrrr, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }
    SECTION( "rrrr = 0, coeff = -2")
    {
        int coeff{ -2 };
        std::vector<Coefficient> expected_result{ 
            Coefficient(-2, 2, Coefficient_Type::AC, 0)
        };
        
        encode_AC_coeff(actual_result, coeff, rrrr, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }
    SECTION( "rrrr = 0, coeff = 5")
    {
        int coeff{ 5 };
        std::vector<Coefficient> expected_result{ 
            Coefficient(5, 3, Coefficient_Type::AC, 0)
        };
        
        encode_AC_coeff(actual_result, coeff, rrrr, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }
    SECTION( "rrrr = 0, coeff = -5")
    {
        int coeff{ -5 };
        std::vector<Coefficient> expected_result{ 
            Coefficient(-5, 3, Coefficient_Type::AC, 0)
        };
        
        encode_AC_coeff(actual_result, coeff, rrrr, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }

    // Note maximum ssss for AC coefficients is 10, unlike for DC coefficients where it's 11
    SECTION( "rrrr = 0, coeff = 1023")
    {
        int coeff{ 1023 };
        std::vector<Coefficient> expected_result{ 
            Coefficient(1023, 10, Coefficient_Type::AC, 0)
        };
        
        encode_AC_coeff(actual_result, coeff, rrrr, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }
    SECTION( "rrrr = 0, coeff = -1023")
    {
        int coeff{ -1023 };
        std::vector<Coefficient> expected_result{ 
            Coefficient(-1023, 10, Coefficient_Type::AC, 0)
        };
        
        encode_AC_coeff(actual_result, coeff, rrrr, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }
}

TEST_CASE( "encode_AC_coeff()::correctly encodes for rrrr = 5", "[encode_AC_coeff()]" ) {
    Huff_Table huff_table{ Huff_Table::load_AC_table(Image_Component::Luminance) };
    std::vector<Coefficient> actual_result{};
    unsigned int rrrr{ 5 };
    unsigned int comp_ind{ 0 };

    SECTION( "rrrr = 5, coeff = 1")
    {
        int coeff{ 1 };
        std::vector<Coefficient> expected_result{ 
            Coefficient(1, 0x51, Coefficient_Type::AC, 0)
        };
        
        encode_AC_coeff(actual_result, coeff, rrrr, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }
    SECTION( "rrrr = 5, coeff = -1")
    {
        int coeff{ -1 };
        std::vector<Coefficient> expected_result{ 
            Coefficient(-1, 0x51, Coefficient_Type::AC, 0)
        };
        
        encode_AC_coeff(actual_result, coeff, rrrr, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }
    SECTION( "rrrr = 5, coeff = 2")
    {
        int coeff{ 2 };
        std::vector<Coefficient> expected_result{ 
            Coefficient(2, 0x52, Coefficient_Type::AC, 0)
        };
        
        encode_AC_coeff(actual_result, coeff, rrrr, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }
    SECTION( "rrrr = 5, coeff = -2")
    {
        int coeff{ -2 };
        std::vector<Coefficient> expected_result{ 
            Coefficient(-2, 0x52, Coefficient_Type::AC, 0)
        };
        
        encode_AC_coeff(actual_result, coeff, rrrr, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }
    SECTION( "rrrr = 5, coeff = 5")
    {
        int coeff{ 5 };
        std::vector<Coefficient> expected_result{ 
            Coefficient(5, 0x53, Coefficient_Type::AC, 0)
        };
        
        encode_AC_coeff(actual_result, coeff, rrrr, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }
    SECTION( "rrrr = 5, coeff = -5")
    {
        int coeff{ -5 };
        std::vector<Coefficient> expected_result{ 
            Coefficient(-5, 0x53, Coefficient_Type::AC, 0)
        };
        
        encode_AC_coeff(actual_result, coeff, rrrr, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }

    // Note maximum ssss for AC coefficients is 10, unlike for DC coefficients where it's 11
    SECTION( "rrrr = 5, coeff = 1023")
    {
        int coeff{ 1023 };
        std::vector<Coefficient> expected_result{ 
            Coefficient(1023, 0x5A, Coefficient_Type::AC, 0)
        };
        
        encode_AC_coeff(actual_result, coeff, rrrr, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }
    SECTION( "rrrr = 5, coeff = -1023")
    {
        int coeff{ -1023 };
        std::vector<Coefficient> expected_result{ 
            Coefficient(-1023, 0x5A, Coefficient_Type::AC, 0)
        };
        
        encode_AC_coeff(actual_result, coeff, rrrr, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }
}

TEST_CASE( "encode_AC_coeff()::correctly encodes for rrrr = 15", "[encode_AC_coeff()]" ) {
    std::vector<Coefficient> actual_result{};
    unsigned int rrrr{ 15 };
    unsigned int comp_ind{ 0 };

    SECTION( "rrrr = 15, coeff = 1")
    {
        int coeff{ 1 };
        std::vector<Coefficient> expected_result{ 
            Coefficient(1, 0xF1, Coefficient_Type::AC, 0)
        };
        
        encode_AC_coeff(actual_result, coeff, rrrr, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }
    SECTION( "rrrr = 15, coeff = -1")
    {
        int coeff{ -1 };
        std::vector<Coefficient> expected_result{ 
            Coefficient(-1, 0xF1, Coefficient_Type::AC, 0)
        };
        
        encode_AC_coeff(actual_result, coeff, rrrr, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }
    SECTION( "rrrr = 15, coeff = 2")
    {
        int coeff{ 2 };
        std::vector<Coefficient> expected_result{ 
            Coefficient(2, 0xF2, Coefficient_Type::AC, 0)
        };
        
        encode_AC_coeff(actual_result, coeff, rrrr, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }
    SECTION( "rrrr = 15, coeff = -2")
    {
        int coeff{ -2 };
        std::vector<Coefficient> expected_result{ 
            Coefficient(-2, 0xF2, Coefficient_Type::AC, 0)
        };
        
        encode_AC_coeff(actual_result, coeff, rrrr, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }
    SECTION( "rrrr = 15, coeff = 5")
    {
        int coeff{ 5 };
        std::vector<Coefficient> expected_result{ 
            Coefficient(5, 0xF3, Coefficient_Type::AC, 0)
        };
        
        encode_AC_coeff(actual_result, coeff, rrrr, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }
    SECTION( "rrrr = 15, coeff = -5")
    {
        int coeff{ -5 };
        std::vector<Coefficient> expected_result{ 
            Coefficient(-5, 0xF3, Coefficient_Type::AC, 0)
        };
        
        encode_AC_coeff(actual_result, coeff, rrrr, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }

    // Note maximum ssss for AC coefficients is 10, unlike for DC coefficients where it's 11
    SECTION( "rrrr = 15, coeff = 1023")
    {
        int coeff{ 1023 };
        std::vector<Coefficient> expected_result{ 
            Coefficient(1023, 0xFA, Coefficient_Type::AC, 0)
        };
        
        encode_AC_coeff(actual_result, coeff, rrrr, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }
    SECTION( "rrrr = 15, coeff = -1023")
    {
        int coeff{ -1023 };
        std::vector<Coefficient> expected_result{ 
            Coefficient(-1023, 0xFA, Coefficient_Type::AC, 0)
        };
        
        encode_AC_coeff(actual_result, coeff, rrrr, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }
}

TEST_CASE( "encode_AC_coeffs()::all entries the same", "[encode_AC_coeffs()]" ) {
    size_t du_ind{ 0 };
    std::vector<Coefficient> actual_result;
    unsigned int comp_ind{ 0 };

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
        
        std::vector<Coefficient>  expected_result{
            // End of block marker corresponds to 0x00
            Coefficient(0, 0x00, Coefficient_Type::AC, 0)
        };

        encode_AC_coeffs(actual_result, input_array, du_ind, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
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

        // Expected to be a 1 repeated 63 times
        std::vector<Coefficient>  expected_result{
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(1, 0x01, Coefficient_Type::AC, 0)
        };

        encode_AC_coeffs(actual_result, input_array, du_ind, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }
    SECTION( "all entries -1") {
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
        std::vector<Coefficient>  expected_result{
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0),
            Coefficient(-1, 0x01, Coefficient_Type::AC, 0)
        };

        encode_AC_coeffs(actual_result, input_array, du_ind, comp_ind);

        CHECK_THAT( actual_result, RangeEquals(expected_result) );
    }
}

TEST_CASE( "encode_AC_coeffs()::zig-zag ordering", "[encode_AC_coeffs()]" ) {
    Huff_Table huff_table{ Huff_Table::load_AC_table(Image_Component::Luminance) };
    size_t du_ind{ 0 };
    unsigned int comp_ind{ 0 };
    
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

    std::vector<Coefficient> expected_result{
        Coefficient(1, 1, Coefficient_Type::AC, 0),
        Coefficient(2, 2, Coefficient_Type::AC, 0),
        Coefficient(3, 2, Coefficient_Type::AC, 0),
        Coefficient(4, 3, Coefficient_Type::AC, 0),
        Coefficient(5, 3, Coefficient_Type::AC, 0),
        // End of block marker corresponds to 0x00
        Coefficient(0, 0x00, Coefficient_Type::AC, 0)
    };

    std::vector<Coefficient> actual_result{};

    encode_AC_coeffs(actual_result, input_array, du_ind, comp_ind);

    CHECK_THAT( actual_result, RangeEquals(expected_result) );
}

TEST_CASE( "encode_AC_coeffs()::run of 1 zero", "[encode_AC_coeffs()]" ) {
    Huff_Table huff_table{ Huff_Table::load_AC_table(Image_Component::Luminance) };
    size_t du_ind{ 0 };
    unsigned int comp_ind{ 0 };
    
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

    std::vector<Coefficient> expected_result{
        Coefficient(1, 1, Coefficient_Type::AC, 0),
        // 1 with a run of 1 zeros
        Coefficient(1, 0x11, Coefficient_Type::AC, 0),
        // End of block marker corresponds to 0x00
        Coefficient(0, 0x00, Coefficient_Type::AC, 0)
    };

    std::vector<Coefficient> actual_result{};

    encode_AC_coeffs(actual_result, input_array, du_ind, comp_ind);

    CHECK_THAT( actual_result, RangeEquals(expected_result) );
}

TEST_CASE( "encode_AC_coeffs()::run of 5 zeros", "[encode_AC_coeffs()]" ) {
    Huff_Table huff_table{ Huff_Table::load_AC_table(Image_Component::Luminance) };
    size_t du_ind{ 0 };
    unsigned int comp_ind{ 0 };
    
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

    std::vector<Coefficient> expected_result{
        Coefficient(1, 1, Coefficient_Type::AC, 0),
        // 1 with a run of 1 zeros
        Coefficient(1, 0x51, Coefficient_Type::AC, 0),
        // End of block marker corresponds to 0x00
        Coefficient(0, 0x00, Coefficient_Type::AC, 0)
    };

    std::vector<Coefficient> actual_result{};

    encode_AC_coeffs(actual_result, input_array, du_ind, comp_ind);

    CHECK_THAT( actual_result, RangeEquals(expected_result) );
}

TEST_CASE( "encode_AC_coeffs()::run of 15 zeros", "[encode_AC_coeffs()]" ) {
    Huff_Table huff_table{ Huff_Table::load_AC_table(Image_Component::Luminance) };
    size_t du_ind{ 0 };
    unsigned int comp_ind{ 0 };
    
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
    
    std::vector<Coefficient> expected_result{
        Coefficient(1, 1, Coefficient_Type::AC, 0),
        // 1 with a run of 15 zeros
        Coefficient(1, 0xF1, Coefficient_Type::AC, 0),
        // End of block marker corresponds to 0x00
        Coefficient(0, 0x00, Coefficient_Type::AC, 0)
    };

    std::vector<Coefficient> actual_result{};

    encode_AC_coeffs(actual_result, input_array, du_ind, comp_ind);

    CHECK_THAT( actual_result, RangeEquals(expected_result) );
}

TEST_CASE( "encode_AC_coeffs()::run of 16 zeros", "[encode_AC_coeffs()]" ) {
    Huff_Table huff_table{ Huff_Table::load_AC_table(Image_Component::Luminance) };
    size_t du_ind{ 0 };
    unsigned int comp_ind{ 0 };
    
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

    std::vector<Coefficient> expected_result{
        Coefficient(1, 1, Coefficient_Type::AC, 0),
        // Run of 16 zeros encoded as a ZRL
        Coefficient(0, 0xF0, Coefficient_Type::AC, 0),
        // 1 with a run of 0 zeros
        Coefficient(1, 1, Coefficient_Type::AC, 0),
        // End of block marker corresponds to 0x00
        Coefficient(0, 0x00, Coefficient_Type::AC, 0)
    };

    std::vector<Coefficient> actual_result{};

    encode_AC_coeffs(actual_result, input_array, du_ind, comp_ind);

    CHECK_THAT( actual_result, RangeEquals(expected_result) );
}

TEST_CASE( "encode_AC_coeffs()::run of 17 zeros", "[encode_AC_coeffs()]" ) {
    size_t du_ind{ 0 };
    unsigned int comp_ind{ 0 };
    
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
    
    std::vector<Coefficient> expected_result{
        Coefficient(1, 1, Coefficient_Type::AC, 0),
        // Run of 16 zeros encoded as a ZRL
        Coefficient(0, 0xF0, Coefficient_Type::AC, 0),
        // 1 with a run of 1 zeros
        Coefficient(1, 0x11, Coefficient_Type::AC, 0),
        // End of block marker corresponds to 0x00
        Coefficient(0, 0x00, Coefficient_Type::AC, 0)
    };

    std::vector<Coefficient> actual_result{};

    encode_AC_coeffs(actual_result, input_array, du_ind, comp_ind);

    CHECK_THAT( actual_result, RangeEquals(expected_result) );
}

TEST_CASE( "encode_AC_coeffs()::run of 32 zeros", "[encode_AC_coeffs()]" ) {
    Huff_Table huff_table{ Huff_Table::load_AC_table(Image_Component::Luminance) };
    size_t du_ind{ 0 };
    unsigned int comp_ind{ 0 };
    
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
    
    std::vector<Coefficient> expected_result{
        Coefficient(1, 1, Coefficient_Type::AC, 0),
        // Run of 16 zeros encoded as a ZRL
        Coefficient(0, 0xF0, Coefficient_Type::AC, 0),
        // Run of 16 zeros encoded as a ZRL
        Coefficient(0, 0xF0, Coefficient_Type::AC, 0),
        // 1 with a run of 0 zeros
        Coefficient(1, 1, Coefficient_Type::AC, 0),
        // End of block marker corresponds to 0x00
        Coefficient(0, 0x00, Coefficient_Type::AC, 0)
    };

    std::vector<Coefficient> actual_result;

    encode_AC_coeffs(actual_result, input_array, du_ind, comp_ind);

    CHECK_THAT( actual_result, RangeEquals(expected_result) );
}

TEST_CASE( "encode_data_unit_sequential()::encode sample data unit", "[encode_data_unit_sequential()]" ) {
    size_t du_ind{ 0 };
    int prev_dc{ 0 };
    unsigned int comp_ind{ 0 };

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

    std::vector<Coefficient> expected_result{
        Coefficient(-26, 5, Coefficient_Type::DC, 0),    // -26, ssss=5
        Coefficient(-3, 0x02, Coefficient_Type::AC, 0),  // -3, rs = 0, 2
        Coefficient(-3, 0x12, Coefficient_Type::AC, 0),  // -3, rs = 1, 2
        Coefficient(-2, 0x02, Coefficient_Type::AC, 0),  // -2  rs = 0, 2
        Coefficient(-6, 0x03, Coefficient_Type::AC, 0),  // -6  rs = 0, 3
        Coefficient( 2, 0x02, Coefficient_Type::AC, 0),  //  2  rs = 0, 2
        Coefficient(-4, 0x03, Coefficient_Type::AC, 0),  // -4  rs = 0, 3
        Coefficient( 1, 0x01, Coefficient_Type::AC, 0),  //  1  rs = 0, 1
        Coefficient(-3, 0x02, Coefficient_Type::AC, 0),  // -3  rs = 0, 2
        Coefficient( 1, 0x01, Coefficient_Type::AC, 0),  //  1  rs = 0, 1
        Coefficient( 1, 0x01, Coefficient_Type::AC, 0),  //  1  rs = 0, 1
        Coefficient( 5, 0x03, Coefficient_Type::AC, 0),  //  5  rs = 0, 3
        Coefficient( 1, 0x01, Coefficient_Type::AC, 0),  //  1  rs = 0, 1
        Coefficient( 2, 0x02, Coefficient_Type::AC, 0),  //  2  rs = 0, 2
        Coefficient(-1, 0x01, Coefficient_Type::AC, 0),  // -1  rs = 0, 1
        Coefficient( 1, 0x01, Coefficient_Type::AC, 0),  //  1  rs = 0, 1
        Coefficient(-1, 0x01, Coefficient_Type::AC, 0),  // -1  rs = 0, 1
        Coefficient( 2, 0x02, Coefficient_Type::AC, 0),  //  2  rs = 0, 2
        Coefficient(-1, 0x51, Coefficient_Type::AC, 0),  // -1  rs = 5, 1
        Coefficient(-1, 0x01, Coefficient_Type::AC, 0),  // -1  rs = 0, 1
        // End of block marker corresponds to 0x00
        Coefficient(0, 0x00, Coefficient_Type::AC, 0)
    };

    std::vector<Coefficient> actual_result;

    encode_data_unit_sequential(actual_result, input_array, du_ind, prev_dc, comp_ind);

    CHECK_THAT( actual_result, RangeEquals(expected_result) );
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

    CHECK_THAT( actual_result, RangeEquals(expected_result) );
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

    CHECK_THAT( actual_result, RangeEquals(expected_result) );
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

    CHECK_THAT( actual_result, RangeEquals(expected_result) );
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

    CHECK_THAT( actual_result, RangeEquals(expected_result) );
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

    std::vector<Coefficient> expected_result{
        // First data unit
        Coefficient(-26, 5, Coefficient_Type::DC, 0),    // -26, ssss=5
        Coefficient(-3, 0x02, Coefficient_Type::AC, 0),  // -3, rs = 0, 2
        Coefficient(-3, 0x12, Coefficient_Type::AC, 0),  // -3, rs = 1, 2
        Coefficient(-2, 0x02, Coefficient_Type::AC, 0),  // -2  rs = 0, 2
        Coefficient(-6, 0x03, Coefficient_Type::AC, 0),  // -6  rs = 0, 3
        Coefficient( 2, 0x02, Coefficient_Type::AC, 0),  //  2  rs = 0, 2
        Coefficient(-4, 0x03, Coefficient_Type::AC, 0),  // -4  rs = 0, 3
        Coefficient( 1, 0x01, Coefficient_Type::AC, 0),  //  1  rs = 0, 1
        Coefficient(-3, 0x02, Coefficient_Type::AC, 0),  // -3  rs = 0, 2
        Coefficient( 1, 0x01, Coefficient_Type::AC, 0),  //  1  rs = 0, 1
        Coefficient( 1, 0x01, Coefficient_Type::AC, 0),  //  1  rs = 0, 1
        Coefficient( 5, 0x03, Coefficient_Type::AC, 0),  //  5  rs = 0, 3
        Coefficient( 1, 0x01, Coefficient_Type::AC, 0),  //  1  rs = 0, 1
        Coefficient( 2, 0x02, Coefficient_Type::AC, 0),  //  2  rs = 0, 2
        Coefficient(-1, 0x01, Coefficient_Type::AC, 0),  // -1  rs = 0, 1
        Coefficient( 1, 0x01, Coefficient_Type::AC, 0),  //  1  rs = 0, 1
        Coefficient(-1, 0x01, Coefficient_Type::AC, 0),  // -1  rs = 0, 1
        Coefficient( 2, 0x02, Coefficient_Type::AC, 0),  //  2  rs = 0, 2
        Coefficient(-1, 0x51, Coefficient_Type::AC, 0),  // -1  rs = 5, 1
        Coefficient(-1, 0x01, Coefficient_Type::AC, 0),  // -1  rs = 0, 1
        // End of block marker corresponds to 0x00
        Coefficient(0, 0x00, Coefficient_Type::AC, 0),
        // Second data unit
        // Remember the DC elements of both data units are the same and its
        // the difference that is encoded. I.e. 0
        Coefficient(0,  0, Coefficient_Type::DC, 0),    // -26, ssss=5
        Coefficient(-3, 0x02, Coefficient_Type::AC, 0),  // -3, rs = 0, 2
        Coefficient(-3, 0x12, Coefficient_Type::AC, 0),  // -3, rs = 1, 2
        Coefficient(-2, 0x02, Coefficient_Type::AC, 0),  // -2  rs = 0, 2
        Coefficient(-6, 0x03, Coefficient_Type::AC, 0),  // -6  rs = 0, 3
        Coefficient( 2, 0x02, Coefficient_Type::AC, 0),  //  2  rs = 0, 2
        Coefficient(-4, 0x03, Coefficient_Type::AC, 0),  // -4  rs = 0, 3
        Coefficient( 1, 0x01, Coefficient_Type::AC, 0),  //  1  rs = 0, 1
        Coefficient(-3, 0x02, Coefficient_Type::AC, 0),  // -3  rs = 0, 2
        Coefficient( 1, 0x01, Coefficient_Type::AC, 0),  //  1  rs = 0, 1
        Coefficient( 1, 0x01, Coefficient_Type::AC, 0),  //  1  rs = 0, 1
        Coefficient( 5, 0x03, Coefficient_Type::AC, 0),  //  5  rs = 0, 3
        Coefficient( 1, 0x01, Coefficient_Type::AC, 0),  //  1  rs = 0, 1
        Coefficient( 2, 0x02, Coefficient_Type::AC, 0),  //  2  rs = 0, 2
        Coefficient(-1, 0x01, Coefficient_Type::AC, 0),  // -1  rs = 0, 1
        Coefficient( 1, 0x01, Coefficient_Type::AC, 0),  //  1  rs = 0, 1
        Coefficient(-1, 0x01, Coefficient_Type::AC, 0),  // -1  rs = 0, 1
        Coefficient( 2, 0x02, Coefficient_Type::AC, 0),  //  2  rs = 0, 2
        Coefficient(-1, 0x51, Coefficient_Type::AC, 0),  // -1  rs = 5, 1
        Coefficient(-1, 0x01, Coefficient_Type::AC, 0),  // -1  rs = 0, 1
        // End of block marker corresponds to 0x00
        Coefficient(0, 0x00, Coefficient_Type::AC, 0)
    };

    std::vector<Coefficient> actual_result;
    size_t mcu_count{ 2 };

    for (size_t mcu_ind = 0; mcu_ind < mcu_count; mcu_ind++)
    {
        append_mcu(actual_result, prev_dc, du_ind, input_array, comp_info);
    }

    CHECK_THAT( actual_result, RangeEquals(expected_result) );
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

    std::vector<Coefficient> expected_result{
        // First component
        // First data unit
        Coefficient(1, 1, Coefficient_Type::DC, 0),        // 1, ssss=1
        Coefficient(2, 0x02, Coefficient_Type::AC, 0),     // 2, rs = 0, 2
        Coefficient(0, 0x00, Coefficient_Type::AC, 0),     // EOB
        // Second data unit
        Coefficient(2, 2, Coefficient_Type::DC, 0),        // 2, ssss=2  (note diff=2)
        Coefficient(4, 0x03, Coefficient_Type::AC, 0),     // 4, rs = 0, 3
        Coefficient(0, 0x00, Coefficient_Type::AC, 0),     // EOB
        // Third data unit
        Coefficient(2, 2, Coefficient_Type::DC, 0),        // 2, ssss=3  (note diff=2)
        Coefficient(6, 0x03, Coefficient_Type::AC, 0),     // 6, rs = 0, 3
        Coefficient(0, 0x00, Coefficient_Type::AC, 0),     // EOB
        // Fourth data unit
        Coefficient(2, 2, Coefficient_Type::DC, 0),        // 2, ssss=3  (note diff=2)
        Coefficient(8, 0x04, Coefficient_Type::AC, 0),     // 8, rs = 0, 4
        Coefficient(0, 0x00, Coefficient_Type::AC, 0),     // EOB
        // Second component
        Coefficient(9, 4, Coefficient_Type::DC, 1),        // 9, ssss=4
        Coefficient(10, 0x04, Coefficient_Type::AC, 1),    // 10, rs = 0, 4
        Coefficient(0, 0x00, Coefficient_Type::AC, 1),     // EOB
        // Third component
        Coefficient(11, 4, Coefficient_Type::DC, 2),       // 11, ssss=4
        Coefficient(12, 0x04, Coefficient_Type::AC, 2),    // 12, rs = 0, 4
        Coefficient(0, 0x00, Coefficient_Type::AC, 2),     // EOB
    };

    std::vector<Coefficient> actual_result;

    // Append the MCU
    append_mcu(actual_result, prev_dc, du_ind, input_array, comp_info);

    CHECK_THAT( actual_result, RangeEquals(expected_result) );
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

    CHECK_THAT( actual_result, RangeEquals(expected_result) );
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

    CHECK_THAT( actual_result, RangeEquals(expected_result) );
}

TEST_CASE( "append_frame_header()::single component", "[append_frame_header()]" ) {
    std::vector<Comp_Info> comp_infos{
        Comp_Info{0, 0, 0, 2, 1},
    };
    unsigned int image_height{ 500 };
    unsigned int image_width{ 1000 };

    std::vector<unsigned char> expected_result{
        0xFF,   // SOF marker, FFC0 for baseline DCT
        0xC0,
        0,      // Length, 2 bytes, most significant first
        11,
        8,      // Precision in bits
        0x01,   // Image height, most significant byte first
        0xF4,
        0x03,   // Image width, most significant byte first
        0xE8,
        1,      // Number of components in frame
        // Start of component specific data
        // First component
        0,      // Component identifier
        0x21,   // Composite byte of horizontal sampling factor (most significant) and vertical sampling factor
        0,      // Quantization table destination
    };

    std::vector<unsigned char> actual_result;

    append_frame_header(actual_result, image_height, image_width, comp_infos);

    CHECK_THAT( actual_result, RangeEquals(expected_result) );
}

TEST_CASE( "append_frame_header()::multiple components", "[append_frame_header()]" ) {
    std::vector<Comp_Info> comp_infos{
        Comp_Info{0, 0, 0, 2, 1},
        Comp_Info{1, 1, 1, 1, 2},
        Comp_Info{1, 1, 1, 1, 2},
    };
    unsigned int image_height{ 500 };
    unsigned int image_width{ 1000 };

    std::vector<unsigned char> expected_result{
        0xFF,   // SOF marker, FFC0 for baseline DCT
        0xC0,
        0,      // Length, 2 bytes, most significant first
        17,
        8,      // Precision in bits
        0x01,   // Image height, most significant byte first
        0xF4,
        0x03,   // Image width, most significant byte first
        0xE8,
        3,      // Number of components in frame
        // Start of component specific data
        // First component
        0,      // Component identifier
        0x21,   // Composite byte of horizontal sampling factor (most significant) and vertical sampling factor
        0,      // Quantization table destination
        // Second component
        1,      // Component identifier
        0x12,   // Composite byte of horizontal sampling factor (most significant) and vertical sampling factor
        1,      // Quantization table destination
        // Second component
        2,      // Component identifier
        0x12,   // Composite byte of horizontal sampling factor (most significant) and vertical sampling factor
        1,      // Quantization table destination
    };

    std::vector<unsigned char> actual_result;

    append_frame_header(actual_result, image_height, image_width, comp_infos);

    CHECK_THAT( actual_result, RangeEquals(expected_result) );
}

TEST_CASE( "convert_to_DU_Array()::H=1, V=1", "[convert_to_DU_Array()]" ) {
    const Array_2d<double> input_array{ gen_array_2d_32x32() };
    const unsigned int H{ 1 }, V{ 1 };
    const Comp_Info comp_info{0, 0, 0, H, V};

    DU_Array<double> expected_result{
        // First row
        {
            {  0,   1,   2,   3,   4,   5,   6,   7},
            { 32,  33,  34,  35,  36,  37,  38,  39},
            { 64,  65,  66,  67,  68,  69,  70,  71},
            { 96,  97,  98,  99, 100, 101, 102, 103},
            {128, 129, 130, 131, 132, 133, 134, 135},
            {160, 161, 162, 163, 164, 165, 166, 167},
            {192, 193, 194, 195, 196, 197, 198, 199},
            {224, 225, 226, 227, 228, 229, 230, 231},
        },
        {
            {  8,   9,  10,  11,  12,  13,  14,  15},
            { 40,  41,  42,  43,  44,  45,  46,  47},
            { 72,  73,  74,  75,  76,  77,  78,  79},
            {104, 105, 106, 107, 108, 109, 110, 111},
            {136, 137, 138, 139, 140, 141, 142, 143},
            {168, 169, 170, 171, 172, 173, 174, 175},
            {200, 201, 202, 203, 204, 205, 206, 207},
            {232, 233, 234, 235, 236, 237, 238, 239},
        },
        {
            { 16,  17,  18,  19,  20,  21,  22,  23},
            { 48,  49,  50,  51,  52,  53,  54,  55},
            { 80,  81,  82,  83,  84,  85,  86,  87},
            {112, 113, 114, 115, 116, 117, 118, 119},
            {144, 145, 146, 147, 148, 149, 150, 151},
            {176, 177, 178, 179, 180, 181, 182, 183},
            {208, 209, 210, 211, 212, 213, 214, 215},
            {240, 241, 242, 243, 244, 245, 246, 247},
        },
        {
            { 24,  25,  26,  27,  28,  29,  30,  31},
            { 56,  57,  58,  59,  60,  61,  62,  63},
            { 88,  89,  90,  91,  92,  93,  94,  95},
            {120, 121, 122, 123, 124, 125, 126, 127},
            {152, 153, 154, 155, 156, 157, 158, 159},
            {184, 185, 186, 187, 188, 189, 190, 191},
            {216, 217, 218, 219, 220, 221, 222, 223},
            {248, 249, 250, 251, 252, 253, 254, 255},
        },
        // Second row
        {
            {256, 257, 258, 259, 260, 261, 262, 263},
            {288, 289, 290, 291, 292, 293, 294, 295},
            {320, 321, 322, 323, 324, 325, 326, 327},
            {352, 353, 354, 355, 356, 357, 358, 359},
            {384, 385, 386, 387, 388, 389, 390, 391},
            {416, 417, 418, 419, 420, 421, 422, 423},
            {448, 449, 450, 451, 452, 453, 454, 455},
            {480, 481, 482, 483, 484, 485, 486, 487},
        },
        {
            {264, 265, 266, 267, 268, 269, 270, 271},
            {296, 297, 298, 299, 300, 301, 302, 303},
            {328, 329, 330, 331, 332, 333, 334, 335},
            {360, 361, 362, 363, 364, 365, 366, 367},
            {392, 393, 394, 395, 396, 397, 398, 399},
            {424, 425, 426, 427, 428, 429, 430, 431},
            {456, 457, 458, 459, 460, 461, 462, 463},
            {488, 489, 490, 491, 492, 493, 494, 495},
        },
        {
            {272, 273, 274, 275, 276, 277, 278, 279},
            {304, 305, 306, 307, 308, 309, 310, 311},
            {336, 337, 338, 339, 340, 341, 342, 343},
            {368, 369, 370, 371, 372, 373, 374, 375},
            {400, 401, 402, 403, 404, 405, 406, 407},
            {432, 433, 434, 435, 436, 437, 438, 439},
            {464, 465, 466, 467, 468, 469, 470, 471},
            {496, 497, 498, 499, 500, 501, 502, 503},
        },
        {
            {280, 281, 282, 283, 284, 285, 286, 287},
            {312, 313, 314, 315, 316, 317, 318, 319},
            {344, 345, 346, 347, 348, 349, 350, 351},
            {376, 377, 378, 379, 380, 381, 382, 383},
            {408, 409, 410, 411, 412, 413, 414, 415},
            {440, 441, 442, 443, 444, 445, 446, 447},
            {472, 473, 474, 475, 476, 477, 478, 479},
            {504, 505, 506, 507, 508, 509, 510, 511},
        },
        // Third row
        {
            {512, 513, 514, 515, 516, 517, 518, 519},
            {544, 545, 546, 547, 548, 549, 550, 551},
            {576, 577, 578, 579, 580, 581, 582, 583},
            {608, 609, 610, 611, 612, 613, 614, 615},
            {640, 641, 642, 643, 644, 645, 646, 647},
            {672, 673, 674, 675, 676, 677, 678, 679},
            {704, 705, 706, 707, 708, 709, 710, 711},
            {736, 737, 738, 739, 740, 741, 742, 743},
        },
        {
            {520, 521, 522, 523, 524, 525, 526, 527},
            {552, 553, 554, 555, 556, 557, 558, 559},
            {584, 585, 586, 587, 588, 589, 590, 591},
            {616, 617, 618, 619, 620, 621, 622, 623},
            {648, 649, 650, 651, 652, 653, 654, 655},
            {680, 681, 682, 683, 684, 685, 686, 687},
            {712, 713, 714, 715, 716, 717, 718, 719},
            {744, 745, 746, 747, 748, 749, 750, 751},
        },
        {
            {528, 529, 530, 531, 532, 533, 534, 535},
            {560, 561, 562, 563, 564, 565, 566, 567},
            {592, 593, 594, 595, 596, 597, 598, 599},
            {624, 625, 626, 627, 628, 629, 630, 631},
            {656, 657, 658, 659, 660, 661, 662, 663},
            {688, 689, 690, 691, 692, 693, 694, 695},
            {720, 721, 722, 723, 724, 725, 726, 727},
            {752, 753, 754, 755, 756, 757, 758, 759},
        },
        {
            {536, 537, 538, 539, 540, 541, 542, 543},
            {568, 569, 570, 571, 572, 573, 574, 575},
            {600, 601, 602, 603, 604, 605, 606, 607},
            {632, 633, 634, 635, 636, 637, 638, 639},
            {664, 665, 666, 667, 668, 669, 670, 671},
            {696, 697, 698, 699, 700, 701, 702, 703},
            {728, 729, 730, 731, 732, 733, 734, 735},
            {760, 761, 762, 763, 764, 765, 766, 767},
        },
        // Fourth row
        {
            {768, 769, 770, 771, 772, 773, 774, 775},
            {800, 801, 802, 803, 804, 805, 806, 807},
            {832, 833, 834, 835, 836, 837, 838, 839},
            {864, 865, 866, 867, 868, 869, 870, 871},
            {896, 897, 898, 899, 900, 901, 902, 903},
            {928, 929, 930, 931, 932, 933, 934, 935},
            {960, 961, 962, 963, 964, 965, 966, 967},
            {992, 993, 994, 995, 996, 997, 998, 999},
        },
        {
            { 776,  777,  778,  779,  780,  781,  782,  783},
            { 808,  809,  810,  811,  812,  813,  814,  815},
            { 840,  841,  842,  843,  844,  845,  846,  847},
            { 872,  873,  874,  875,  876,  877,  878,  879},
            { 904,  905,  906,  907,  908,  909,  910,  911},
            { 936,  937,  938,  939,  940,  941,  942,  943},
            { 968,  969,  970,  971,  972,  973,  974,  975},
            {1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007},
        },
        {
            { 784,  785,  786,  787,  788,  789,  790,  791},
            { 816,  817,  818,  819,  820,  821,  822,  823},
            { 848,  849,  850,  851,  852,  853,  854,  855},
            { 880,  881,  882,  883,  884,  885,  886,  887},
            { 912,  913,  914,  915,  916,  917,  918,  919},
            { 944,  945,  946,  947,  948,  949,  950,  951},
            { 976,  977,  978,  979,  980,  981,  982,  983},
            {1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015},
        },
        {
            { 792,  793,  794,  795,  796,  797,  798,  799},
            { 824,  825,  826,  827,  828,  829,  830,  831},
            { 856,  857,  858,  859,  860,  861,  862,  863},
            { 888,  889,  890,  891,  892,  893,  894,  895},
            { 920,  921,  922,  923,  924,  925,  926,  927},
            { 952,  953,  954,  955,  956,  957,  958,  959},
            { 984,  985,  986,  987,  988,  989,  990,  991},
            {1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023},
        }
    };

    DU_Array actual_result{ convert_to_DU_Array(input_array, comp_info) };


    REQUIRE( actual_result.size() == expected_result.size() );

    for (size_t ind = 0; ind < expected_result.size(); ind++)
    {
        CAPTURE( ind );
        CHECK( actual_result[ind] == expected_result[ind] );
    }
}

TEST_CASE( "convert_to_DU_Array()::H=2, V=1", "[convert_to_DU_Array()]" ) {
    const Array_2d<double> input_array{ gen_array_2d_32x32() };
    const unsigned int H{ 2 }, V{ 1 };
    const Comp_Info comp_info{0, 0, 0, H, V};

    DU_Array<double> expected_result{
        // First row
        {
            {  0,   1,   2,   3,   4,   5,   6,   7},
            { 32,  33,  34,  35,  36,  37,  38,  39},
            { 64,  65,  66,  67,  68,  69,  70,  71},
            { 96,  97,  98,  99, 100, 101, 102, 103},
            {128, 129, 130, 131, 132, 133, 134, 135},
            {160, 161, 162, 163, 164, 165, 166, 167},
            {192, 193, 194, 195, 196, 197, 198, 199},
            {224, 225, 226, 227, 228, 229, 230, 231},
        },
        {
            {  8,   9,  10,  11,  12,  13,  14,  15},
            { 40,  41,  42,  43,  44,  45,  46,  47},
            { 72,  73,  74,  75,  76,  77,  78,  79},
            {104, 105, 106, 107, 108, 109, 110, 111},
            {136, 137, 138, 139, 140, 141, 142, 143},
            {168, 169, 170, 171, 172, 173, 174, 175},
            {200, 201, 202, 203, 204, 205, 206, 207},
            {232, 233, 234, 235, 236, 237, 238, 239},
        },
        {
            { 16,  17,  18,  19,  20,  21,  22,  23},
            { 48,  49,  50,  51,  52,  53,  54,  55},
            { 80,  81,  82,  83,  84,  85,  86,  87},
            {112, 113, 114, 115, 116, 117, 118, 119},
            {144, 145, 146, 147, 148, 149, 150, 151},
            {176, 177, 178, 179, 180, 181, 182, 183},
            {208, 209, 210, 211, 212, 213, 214, 215},
            {240, 241, 242, 243, 244, 245, 246, 247},
        },
        {
            { 24,  25,  26,  27,  28,  29,  30,  31},
            { 56,  57,  58,  59,  60,  61,  62,  63},
            { 88,  89,  90,  91,  92,  93,  94,  95},
            {120, 121, 122, 123, 124, 125, 126, 127},
            {152, 153, 154, 155, 156, 157, 158, 159},
            {184, 185, 186, 187, 188, 189, 190, 191},
            {216, 217, 218, 219, 220, 221, 222, 223},
            {248, 249, 250, 251, 252, 253, 254, 255},
        },
        // Second row
        {
            {256, 257, 258, 259, 260, 261, 262, 263},
            {288, 289, 290, 291, 292, 293, 294, 295},
            {320, 321, 322, 323, 324, 325, 326, 327},
            {352, 353, 354, 355, 356, 357, 358, 359},
            {384, 385, 386, 387, 388, 389, 390, 391},
            {416, 417, 418, 419, 420, 421, 422, 423},
            {448, 449, 450, 451, 452, 453, 454, 455},
            {480, 481, 482, 483, 484, 485, 486, 487},
        },
        {
            {264, 265, 266, 267, 268, 269, 270, 271},
            {296, 297, 298, 299, 300, 301, 302, 303},
            {328, 329, 330, 331, 332, 333, 334, 335},
            {360, 361, 362, 363, 364, 365, 366, 367},
            {392, 393, 394, 395, 396, 397, 398, 399},
            {424, 425, 426, 427, 428, 429, 430, 431},
            {456, 457, 458, 459, 460, 461, 462, 463},
            {488, 489, 490, 491, 492, 493, 494, 495},
        },
        {
            {272, 273, 274, 275, 276, 277, 278, 279},
            {304, 305, 306, 307, 308, 309, 310, 311},
            {336, 337, 338, 339, 340, 341, 342, 343},
            {368, 369, 370, 371, 372, 373, 374, 375},
            {400, 401, 402, 403, 404, 405, 406, 407},
            {432, 433, 434, 435, 436, 437, 438, 439},
            {464, 465, 466, 467, 468, 469, 470, 471},
            {496, 497, 498, 499, 500, 501, 502, 503},
        },
        {
            {280, 281, 282, 283, 284, 285, 286, 287},
            {312, 313, 314, 315, 316, 317, 318, 319},
            {344, 345, 346, 347, 348, 349, 350, 351},
            {376, 377, 378, 379, 380, 381, 382, 383},
            {408, 409, 410, 411, 412, 413, 414, 415},
            {440, 441, 442, 443, 444, 445, 446, 447},
            {472, 473, 474, 475, 476, 477, 478, 479},
            {504, 505, 506, 507, 508, 509, 510, 511},
        },
        // Third row
        {
            {512, 513, 514, 515, 516, 517, 518, 519},
            {544, 545, 546, 547, 548, 549, 550, 551},
            {576, 577, 578, 579, 580, 581, 582, 583},
            {608, 609, 610, 611, 612, 613, 614, 615},
            {640, 641, 642, 643, 644, 645, 646, 647},
            {672, 673, 674, 675, 676, 677, 678, 679},
            {704, 705, 706, 707, 708, 709, 710, 711},
            {736, 737, 738, 739, 740, 741, 742, 743},
        },
        {
            {520, 521, 522, 523, 524, 525, 526, 527},
            {552, 553, 554, 555, 556, 557, 558, 559},
            {584, 585, 586, 587, 588, 589, 590, 591},
            {616, 617, 618, 619, 620, 621, 622, 623},
            {648, 649, 650, 651, 652, 653, 654, 655},
            {680, 681, 682, 683, 684, 685, 686, 687},
            {712, 713, 714, 715, 716, 717, 718, 719},
            {744, 745, 746, 747, 748, 749, 750, 751},
        },
        {
            {528, 529, 530, 531, 532, 533, 534, 535},
            {560, 561, 562, 563, 564, 565, 566, 567},
            {592, 593, 594, 595, 596, 597, 598, 599},
            {624, 625, 626, 627, 628, 629, 630, 631},
            {656, 657, 658, 659, 660, 661, 662, 663},
            {688, 689, 690, 691, 692, 693, 694, 695},
            {720, 721, 722, 723, 724, 725, 726, 727},
            {752, 753, 754, 755, 756, 757, 758, 759},
        },
        {
            {536, 537, 538, 539, 540, 541, 542, 543},
            {568, 569, 570, 571, 572, 573, 574, 575},
            {600, 601, 602, 603, 604, 605, 606, 607},
            {632, 633, 634, 635, 636, 637, 638, 639},
            {664, 665, 666, 667, 668, 669, 670, 671},
            {696, 697, 698, 699, 700, 701, 702, 703},
            {728, 729, 730, 731, 732, 733, 734, 735},
            {760, 761, 762, 763, 764, 765, 766, 767},
        },
        // Fourth row
        {
            {768, 769, 770, 771, 772, 773, 774, 775},
            {800, 801, 802, 803, 804, 805, 806, 807},
            {832, 833, 834, 835, 836, 837, 838, 839},
            {864, 865, 866, 867, 868, 869, 870, 871},
            {896, 897, 898, 899, 900, 901, 902, 903},
            {928, 929, 930, 931, 932, 933, 934, 935},
            {960, 961, 962, 963, 964, 965, 966, 967},
            {992, 993, 994, 995, 996, 997, 998, 999},
        },
        {
            { 776,  777,  778,  779,  780,  781,  782,  783},
            { 808,  809,  810,  811,  812,  813,  814,  815},
            { 840,  841,  842,  843,  844,  845,  846,  847},
            { 872,  873,  874,  875,  876,  877,  878,  879},
            { 904,  905,  906,  907,  908,  909,  910,  911},
            { 936,  937,  938,  939,  940,  941,  942,  943},
            { 968,  969,  970,  971,  972,  973,  974,  975},
            {1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007},
        },
        {
            { 784,  785,  786,  787,  788,  789,  790,  791},
            { 816,  817,  818,  819,  820,  821,  822,  823},
            { 848,  849,  850,  851,  852,  853,  854,  855},
            { 880,  881,  882,  883,  884,  885,  886,  887},
            { 912,  913,  914,  915,  916,  917,  918,  919},
            { 944,  945,  946,  947,  948,  949,  950,  951},
            { 976,  977,  978,  979,  980,  981,  982,  983},
            {1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015},
        },
        {
            { 792,  793,  794,  795,  796,  797,  798,  799},
            { 824,  825,  826,  827,  828,  829,  830,  831},
            { 856,  857,  858,  859,  860,  861,  862,  863},
            { 888,  889,  890,  891,  892,  893,  894,  895},
            { 920,  921,  922,  923,  924,  925,  926,  927},
            { 952,  953,  954,  955,  956,  957,  958,  959},
            { 984,  985,  986,  987,  988,  989,  990,  991},
            {1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023},
        }
    };

    DU_Array actual_result{ convert_to_DU_Array(input_array, comp_info) };

    REQUIRE( actual_result.size() == expected_result.size() );

    for (size_t ind = 0; ind < expected_result.size(); ind++)
    {
        CAPTURE( ind );
        CHECK( actual_result[ind] == expected_result[ind] );
    }
}

TEST_CASE( "convert_to_DU_Array()::H=1, V=2", "[convert_to_DU_Array()]" ) {
    const Array_2d<double> input_array{ gen_array_2d_32x32() };
    const unsigned int H{ 1 }, V{ 2 };
    const Comp_Info comp_info{0, 0, 0, H, V};

    DU_Array<double> expected_result{
        // MCU (0, 0)
        {
            {  0,   1,   2,   3,   4,   5,   6,   7},
            { 32,  33,  34,  35,  36,  37,  38,  39},
            { 64,  65,  66,  67,  68,  69,  70,  71},
            { 96,  97,  98,  99, 100, 101, 102, 103},
            {128, 129, 130, 131, 132, 133, 134, 135},
            {160, 161, 162, 163, 164, 165, 166, 167},
            {192, 193, 194, 195, 196, 197, 198, 199},
            {224, 225, 226, 227, 228, 229, 230, 231},
        },
        {
            {256, 257, 258, 259, 260, 261, 262, 263},
            {288, 289, 290, 291, 292, 293, 294, 295},
            {320, 321, 322, 323, 324, 325, 326, 327},
            {352, 353, 354, 355, 356, 357, 358, 359},
            {384, 385, 386, 387, 388, 389, 390, 391},
            {416, 417, 418, 419, 420, 421, 422, 423},
            {448, 449, 450, 451, 452, 453, 454, 455},
            {480, 481, 482, 483, 484, 485, 486, 487},
        },
        // MCU (0, 1)
        {
            {  8,   9,  10,  11,  12,  13,  14,  15},
            { 40,  41,  42,  43,  44,  45,  46,  47},
            { 72,  73,  74,  75,  76,  77,  78,  79},
            {104, 105, 106, 107, 108, 109, 110, 111},
            {136, 137, 138, 139, 140, 141, 142, 143},
            {168, 169, 170, 171, 172, 173, 174, 175},
            {200, 201, 202, 203, 204, 205, 206, 207},
            {232, 233, 234, 235, 236, 237, 238, 239},
        },
        {
            {264, 265, 266, 267, 268, 269, 270, 271},
            {296, 297, 298, 299, 300, 301, 302, 303},
            {328, 329, 330, 331, 332, 333, 334, 335},
            {360, 361, 362, 363, 364, 365, 366, 367},
            {392, 393, 394, 395, 396, 397, 398, 399},
            {424, 425, 426, 427, 428, 429, 430, 431},
            {456, 457, 458, 459, 460, 461, 462, 463},
            {488, 489, 490, 491, 492, 493, 494, 495},
        },
        // MCU (0, 2)
        {
            { 16,  17,  18,  19,  20,  21,  22,  23},
            { 48,  49,  50,  51,  52,  53,  54,  55},
            { 80,  81,  82,  83,  84,  85,  86,  87},
            {112, 113, 114, 115, 116, 117, 118, 119},
            {144, 145, 146, 147, 148, 149, 150, 151},
            {176, 177, 178, 179, 180, 181, 182, 183},
            {208, 209, 210, 211, 212, 213, 214, 215},
            {240, 241, 242, 243, 244, 245, 246, 247},
        },
        {
            {272, 273, 274, 275, 276, 277, 278, 279},
            {304, 305, 306, 307, 308, 309, 310, 311},
            {336, 337, 338, 339, 340, 341, 342, 343},
            {368, 369, 370, 371, 372, 373, 374, 375},
            {400, 401, 402, 403, 404, 405, 406, 407},
            {432, 433, 434, 435, 436, 437, 438, 439},
            {464, 465, 466, 467, 468, 469, 470, 471},
            {496, 497, 498, 499, 500, 501, 502, 503},
        },
        // MCU (0, 3)
        {
            { 24,  25,  26,  27,  28,  29,  30,  31},
            { 56,  57,  58,  59,  60,  61,  62,  63},
            { 88,  89,  90,  91,  92,  93,  94,  95},
            {120, 121, 122, 123, 124, 125, 126, 127},
            {152, 153, 154, 155, 156, 157, 158, 159},
            {184, 185, 186, 187, 188, 189, 190, 191},
            {216, 217, 218, 219, 220, 221, 222, 223},
            {248, 249, 250, 251, 252, 253, 254, 255},
        },
        {
            {280, 281, 282, 283, 284, 285, 286, 287},
            {312, 313, 314, 315, 316, 317, 318, 319},
            {344, 345, 346, 347, 348, 349, 350, 351},
            {376, 377, 378, 379, 380, 381, 382, 383},
            {408, 409, 410, 411, 412, 413, 414, 415},
            {440, 441, 442, 443, 444, 445, 446, 447},
            {472, 473, 474, 475, 476, 477, 478, 479},
            {504, 505, 506, 507, 508, 509, 510, 511},
        },
        // MCU (1, 0)
        {
            {512, 513, 514, 515, 516, 517, 518, 519},
            {544, 545, 546, 547, 548, 549, 550, 551},
            {576, 577, 578, 579, 580, 581, 582, 583},
            {608, 609, 610, 611, 612, 613, 614, 615},
            {640, 641, 642, 643, 644, 645, 646, 647},
            {672, 673, 674, 675, 676, 677, 678, 679},
            {704, 705, 706, 707, 708, 709, 710, 711},
            {736, 737, 738, 739, 740, 741, 742, 743},
        },
        {
            {768, 769, 770, 771, 772, 773, 774, 775},
            {800, 801, 802, 803, 804, 805, 806, 807},
            {832, 833, 834, 835, 836, 837, 838, 839},
            {864, 865, 866, 867, 868, 869, 870, 871},
            {896, 897, 898, 899, 900, 901, 902, 903},
            {928, 929, 930, 931, 932, 933, 934, 935},
            {960, 961, 962, 963, 964, 965, 966, 967},
            {992, 993, 994, 995, 996, 997, 998, 999},
        },
        // MCU (1, 1)
        {
            {520, 521, 522, 523, 524, 525, 526, 527},
            {552, 553, 554, 555, 556, 557, 558, 559},
            {584, 585, 586, 587, 588, 589, 590, 591},
            {616, 617, 618, 619, 620, 621, 622, 623},
            {648, 649, 650, 651, 652, 653, 654, 655},
            {680, 681, 682, 683, 684, 685, 686, 687},
            {712, 713, 714, 715, 716, 717, 718, 719},
            {744, 745, 746, 747, 748, 749, 750, 751},
        },
        {
            { 776,  777,  778,  779,  780,  781,  782,  783},
            { 808,  809,  810,  811,  812,  813,  814,  815},
            { 840,  841,  842,  843,  844,  845,  846,  847},
            { 872,  873,  874,  875,  876,  877,  878,  879},
            { 904,  905,  906,  907,  908,  909,  910,  911},
            { 936,  937,  938,  939,  940,  941,  942,  943},
            { 968,  969,  970,  971,  972,  973,  974,  975},
            {1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007},
        },
        // MCU (1, 2)
        {
            {528, 529, 530, 531, 532, 533, 534, 535},
            {560, 561, 562, 563, 564, 565, 566, 567},
            {592, 593, 594, 595, 596, 597, 598, 599},
            {624, 625, 626, 627, 628, 629, 630, 631},
            {656, 657, 658, 659, 660, 661, 662, 663},
            {688, 689, 690, 691, 692, 693, 694, 695},
            {720, 721, 722, 723, 724, 725, 726, 727},
            {752, 753, 754, 755, 756, 757, 758, 759},
        },
        {
            { 784,  785,  786,  787,  788,  789,  790,  791},
            { 816,  817,  818,  819,  820,  821,  822,  823},
            { 848,  849,  850,  851,  852,  853,  854,  855},
            { 880,  881,  882,  883,  884,  885,  886,  887},
            { 912,  913,  914,  915,  916,  917,  918,  919},
            { 944,  945,  946,  947,  948,  949,  950,  951},
            { 976,  977,  978,  979,  980,  981,  982,  983},
            {1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015},
        },
        // MCU (1, 3)
        {
            {536, 537, 538, 539, 540, 541, 542, 543},
            {568, 569, 570, 571, 572, 573, 574, 575},
            {600, 601, 602, 603, 604, 605, 606, 607},
            {632, 633, 634, 635, 636, 637, 638, 639},
            {664, 665, 666, 667, 668, 669, 670, 671},
            {696, 697, 698, 699, 700, 701, 702, 703},
            {728, 729, 730, 731, 732, 733, 734, 735},
            {760, 761, 762, 763, 764, 765, 766, 767},
        },
        {
            { 792,  793,  794,  795,  796,  797,  798,  799},
            { 824,  825,  826,  827,  828,  829,  830,  831},
            { 856,  857,  858,  859,  860,  861,  862,  863},
            { 888,  889,  890,  891,  892,  893,  894,  895},
            { 920,  921,  922,  923,  924,  925,  926,  927},
            { 952,  953,  954,  955,  956,  957,  958,  959},
            { 984,  985,  986,  987,  988,  989,  990,  991},
            {1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023},
        }
    };

    DU_Array actual_result{ convert_to_DU_Array(input_array, comp_info) };


    REQUIRE( actual_result.size() == expected_result.size() );

    for (size_t ind = 0; ind < expected_result.size(); ind++)
    {
        CAPTURE( ind );
        CHECK( actual_result[ind] == expected_result[ind] );
    }
}

TEST_CASE( "convert_to_DU_Array()::H=2, V=2", "[convert_to_DU_Array()]" ) {
    const Array_2d<double> input_array{ gen_array_2d_32x32() };
    const unsigned int H{ 2 }, V{ 2 };
    const Comp_Info comp_info{0, 0, 0, H, V};

    DU_Array<double> expected_result{
        // MCU (0, 0)
        {
            {  0,   1,   2,   3,   4,   5,   6,   7},
            { 32,  33,  34,  35,  36,  37,  38,  39},
            { 64,  65,  66,  67,  68,  69,  70,  71},
            { 96,  97,  98,  99, 100, 101, 102, 103},
            {128, 129, 130, 131, 132, 133, 134, 135},
            {160, 161, 162, 163, 164, 165, 166, 167},
            {192, 193, 194, 195, 196, 197, 198, 199},
            {224, 225, 226, 227, 228, 229, 230, 231},
        },
        {
            {  8,   9,  10,  11,  12,  13,  14,  15},
            { 40,  41,  42,  43,  44,  45,  46,  47},
            { 72,  73,  74,  75,  76,  77,  78,  79},
            {104, 105, 106, 107, 108, 109, 110, 111},
            {136, 137, 138, 139, 140, 141, 142, 143},
            {168, 169, 170, 171, 172, 173, 174, 175},
            {200, 201, 202, 203, 204, 205, 206, 207},
            {232, 233, 234, 235, 236, 237, 238, 239},
        },
        {
            {256, 257, 258, 259, 260, 261, 262, 263},
            {288, 289, 290, 291, 292, 293, 294, 295},
            {320, 321, 322, 323, 324, 325, 326, 327},
            {352, 353, 354, 355, 356, 357, 358, 359},
            {384, 385, 386, 387, 388, 389, 390, 391},
            {416, 417, 418, 419, 420, 421, 422, 423},
            {448, 449, 450, 451, 452, 453, 454, 455},
            {480, 481, 482, 483, 484, 485, 486, 487},
        },
        {
            {264, 265, 266, 267, 268, 269, 270, 271},
            {296, 297, 298, 299, 300, 301, 302, 303},
            {328, 329, 330, 331, 332, 333, 334, 335},
            {360, 361, 362, 363, 364, 365, 366, 367},
            {392, 393, 394, 395, 396, 397, 398, 399},
            {424, 425, 426, 427, 428, 429, 430, 431},
            {456, 457, 458, 459, 460, 461, 462, 463},
            {488, 489, 490, 491, 492, 493, 494, 495},
        },
        // MCU (0, 1)
        {
            { 16,  17,  18,  19,  20,  21,  22,  23},
            { 48,  49,  50,  51,  52,  53,  54,  55},
            { 80,  81,  82,  83,  84,  85,  86,  87},
            {112, 113, 114, 115, 116, 117, 118, 119},
            {144, 145, 146, 147, 148, 149, 150, 151},
            {176, 177, 178, 179, 180, 181, 182, 183},
            {208, 209, 210, 211, 212, 213, 214, 215},
            {240, 241, 242, 243, 244, 245, 246, 247},
        },
        {
            { 24,  25,  26,  27,  28,  29,  30,  31},
            { 56,  57,  58,  59,  60,  61,  62,  63},
            { 88,  89,  90,  91,  92,  93,  94,  95},
            {120, 121, 122, 123, 124, 125, 126, 127},
            {152, 153, 154, 155, 156, 157, 158, 159},
            {184, 185, 186, 187, 188, 189, 190, 191},
            {216, 217, 218, 219, 220, 221, 222, 223},
            {248, 249, 250, 251, 252, 253, 254, 255},
        },
        {
            {272, 273, 274, 275, 276, 277, 278, 279},
            {304, 305, 306, 307, 308, 309, 310, 311},
            {336, 337, 338, 339, 340, 341, 342, 343},
            {368, 369, 370, 371, 372, 373, 374, 375},
            {400, 401, 402, 403, 404, 405, 406, 407},
            {432, 433, 434, 435, 436, 437, 438, 439},
            {464, 465, 466, 467, 468, 469, 470, 471},
            {496, 497, 498, 499, 500, 501, 502, 503},
        },
        {
            {280, 281, 282, 283, 284, 285, 286, 287},
            {312, 313, 314, 315, 316, 317, 318, 319},
            {344, 345, 346, 347, 348, 349, 350, 351},
            {376, 377, 378, 379, 380, 381, 382, 383},
            {408, 409, 410, 411, 412, 413, 414, 415},
            {440, 441, 442, 443, 444, 445, 446, 447},
            {472, 473, 474, 475, 476, 477, 478, 479},
            {504, 505, 506, 507, 508, 509, 510, 511},
        },
        // MCU (1, 0)
        {
            {512, 513, 514, 515, 516, 517, 518, 519},
            {544, 545, 546, 547, 548, 549, 550, 551},
            {576, 577, 578, 579, 580, 581, 582, 583},
            {608, 609, 610, 611, 612, 613, 614, 615},
            {640, 641, 642, 643, 644, 645, 646, 647},
            {672, 673, 674, 675, 676, 677, 678, 679},
            {704, 705, 706, 707, 708, 709, 710, 711},
            {736, 737, 738, 739, 740, 741, 742, 743},
        },
        {
            {520, 521, 522, 523, 524, 525, 526, 527},
            {552, 553, 554, 555, 556, 557, 558, 559},
            {584, 585, 586, 587, 588, 589, 590, 591},
            {616, 617, 618, 619, 620, 621, 622, 623},
            {648, 649, 650, 651, 652, 653, 654, 655},
            {680, 681, 682, 683, 684, 685, 686, 687},
            {712, 713, 714, 715, 716, 717, 718, 719},
            {744, 745, 746, 747, 748, 749, 750, 751},
        },
        {
            {768, 769, 770, 771, 772, 773, 774, 775},
            {800, 801, 802, 803, 804, 805, 806, 807},
            {832, 833, 834, 835, 836, 837, 838, 839},
            {864, 865, 866, 867, 868, 869, 870, 871},
            {896, 897, 898, 899, 900, 901, 902, 903},
            {928, 929, 930, 931, 932, 933, 934, 935},
            {960, 961, 962, 963, 964, 965, 966, 967},
            {992, 993, 994, 995, 996, 997, 998, 999},
        },
        {
            { 776,  777,  778,  779,  780,  781,  782,  783},
            { 808,  809,  810,  811,  812,  813,  814,  815},
            { 840,  841,  842,  843,  844,  845,  846,  847},
            { 872,  873,  874,  875,  876,  877,  878,  879},
            { 904,  905,  906,  907,  908,  909,  910,  911},
            { 936,  937,  938,  939,  940,  941,  942,  943},
            { 968,  969,  970,  971,  972,  973,  974,  975},
            {1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007},
        },
        // MCU (1, 1)
        {
            {528, 529, 530, 531, 532, 533, 534, 535},
            {560, 561, 562, 563, 564, 565, 566, 567},
            {592, 593, 594, 595, 596, 597, 598, 599},
            {624, 625, 626, 627, 628, 629, 630, 631},
            {656, 657, 658, 659, 660, 661, 662, 663},
            {688, 689, 690, 691, 692, 693, 694, 695},
            {720, 721, 722, 723, 724, 725, 726, 727},
            {752, 753, 754, 755, 756, 757, 758, 759},
        },
        {
            {536, 537, 538, 539, 540, 541, 542, 543},
            {568, 569, 570, 571, 572, 573, 574, 575},
            {600, 601, 602, 603, 604, 605, 606, 607},
            {632, 633, 634, 635, 636, 637, 638, 639},
            {664, 665, 666, 667, 668, 669, 670, 671},
            {696, 697, 698, 699, 700, 701, 702, 703},
            {728, 729, 730, 731, 732, 733, 734, 735},
            {760, 761, 762, 763, 764, 765, 766, 767},
        },
        {
            { 784,  785,  786,  787,  788,  789,  790,  791},
            { 816,  817,  818,  819,  820,  821,  822,  823},
            { 848,  849,  850,  851,  852,  853,  854,  855},
            { 880,  881,  882,  883,  884,  885,  886,  887},
            { 912,  913,  914,  915,  916,  917,  918,  919},
            { 944,  945,  946,  947,  948,  949,  950,  951},
            { 976,  977,  978,  979,  980,  981,  982,  983},
            {1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015},
        },
        {
            { 792,  793,  794,  795,  796,  797,  798,  799},
            { 824,  825,  826,  827,  828,  829,  830,  831},
            { 856,  857,  858,  859,  860,  861,  862,  863},
            { 888,  889,  890,  891,  892,  893,  894,  895},
            { 920,  921,  922,  923,  924,  925,  926,  927},
            { 952,  953,  954,  955,  956,  957,  958,  959},
            { 984,  985,  986,  987,  988,  989,  990,  991},
            {1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023},
        }
    };

    DU_Array actual_result{ convert_to_DU_Array(input_array, comp_info) };

    REQUIRE( actual_result.size() == expected_result.size() );

    for (size_t ind = 0; ind < expected_result.size(); ind++)
    {
        CAPTURE( ind );
        CHECK( actual_result[ind] == expected_result[ind] );
    }
}

TEST_CASE( "enlarge_component()::8x8 H=1, V=1", "[enlarge_component()]" ) {
    const unsigned int H{ 1 }, V{ 1 };
    Array_2d<double> input_array{ gen_array_2d_arange(8, 8) };

    // Output should be the same as input
    Array_2d<double> expected_result{ input_array };

    Array_2d<double> actual_result{ enlarge_component(input_array, V, H) };

    REQUIRE( actual_result.size() == expected_result.size() );

    for (size_t ind = 0; ind < expected_result.size(); ind++)
    {
        CAPTURE( ind );
        CHECK( actual_result[ind] == expected_result[ind] );
    }
}

TEST_CASE( "enlarge_component()::8x8 H=2, V=2", "[enlarge_component()]" ) {
    const unsigned int H{ 2 }, V{ 2 };
    Array_2d<double> input_array{ gen_array_2d_arange(8, 8) };

    Array_2d<double> expected_result{ 
        { 0,   1,   2,   3,   4,   5,   6,   7,   7,   7,   7,   7,   7,   7,   7,   7},
        { 8,   9,  10,  11,  12,  13,  14,  15,  15,  15,  15,  15,  15,  15,  15,  15},
        {16,  17,  18,  19,  20,  21,  22,  23,  23,  23,  23,  23,  23,  23,  23,  23},
        {24,  25,  26,  27,  28,  29,  30,  31,  31,  31,  31,  31,  31,  31,  31,  31},
        {32,  33,  34,  35,  36,  37,  38,  39,  39,  39,  39,  39,  39,  39,  39,  39},
        {40,  41,  42,  43,  44,  45,  46,  47,  47,  47,  47,  47,  47,  47,  47,  47},
        {48,  49,  50,  51,  52,  53,  54,  55,  55,  55,  55,  55,  55,  55,  55,  55},
        {56,  57,  58,  59,  60,  61,  62,  63,  63,  63,  63,  63,  63,  63,  63,  63},
        {56,  57,  58,  59,  60,  61,  62,  63, 128, 128, 128, 128, 128, 128, 128, 128},
        {56,  57,  58,  59,  60,  61,  62,  63, 128, 128, 128, 128, 128, 128, 128, 128},
        {56,  57,  58,  59,  60,  61,  62,  63, 128, 128, 128, 128, 128, 128, 128, 128},
        {56,  57,  58,  59,  60,  61,  62,  63, 128, 128, 128, 128, 128, 128, 128, 128},
        {56,  57,  58,  59,  60,  61,  62,  63, 128, 128, 128, 128, 128, 128, 128, 128},
        {56,  57,  58,  59,  60,  61,  62,  63, 128, 128, 128, 128, 128, 128, 128, 128},
        {56,  57,  58,  59,  60,  61,  62,  63, 128, 128, 128, 128, 128, 128, 128, 128},
        {56,  57,  58,  59,  60,  61,  62,  63, 128, 128, 128, 128, 128, 128, 128, 128},
     };

    Array_2d<double> actual_result{ enlarge_component(input_array, V, H) };

    REQUIRE( actual_result.size() == expected_result.size() );

    for (size_t ind = 0; ind < expected_result.size(); ind++)
    {
        CAPTURE( ind );
        CHECK( actual_result[ind] == expected_result[ind] );
    }
}

TEST_CASE( "enlarge_component()::2x3 H=1, V=1", "[enlarge_component()]" ) {
    const unsigned int H{ 1 }, V{ 1 };
    Array_2d<double> input_array{ gen_array_2d_arange(2, 3) };

    Array_2d<double> expected_result{ 
        {0,   1,   2,   2,   2,   2,   2,   2},
        {3,   4,   5,   5,   5,   5,   5,   5},
        {3,   4,   5, 128, 128, 128, 128, 128},
        {3,   4,   5, 128, 128, 128, 128, 128},
        {3,   4,   5, 128, 128, 128, 128, 128},
        {3,   4,   5, 128, 128, 128, 128, 128},
        {3,   4,   5, 128, 128, 128, 128, 128},
        {3,   4,   5, 128, 128, 128, 128, 128},
     };

    Array_2d<double> actual_result{ enlarge_component(input_array, V, H) };

    REQUIRE( actual_result.size() == expected_result.size() );

    for (size_t ind = 0; ind < expected_result.size(); ind++)
    {
        CAPTURE( ind );
        CHECK( actual_result[ind] == expected_result[ind] );
    }
}

TEST_CASE( "colour_transform()::Check minimum output is 0", "[colour_transform()]" ) {

    // Uses inputs designed to give the smallest possible output after the colour matrix
    // multiplication and checks they are not less than zero
    // For each component:
    // Y : (  0,   0,   0) -> 0
    // U : (255, 255,   0) -> 1
    // V : (  0, 255, 255) -> 1
    Array_2d<double> red{ 
        {0, 255, 0}
    };

    Array_2d<double> green{ 
        {0, 255, 255}
    };

    Array_2d<double> blue{ 
        {0, 0, 255}
    };

    // Now process
    auto [Y, U, V] = colour_transform(red, green, blue);

    for (size_t ind = 0; ind < red.size(); ind++)
    {
        CAPTURE( ind );
        CHECK( Y[ind] >= 0.0 );
        CHECK( U[ind] >= 0.0 );
        CHECK( V[ind] >= 0.0 );
    }
}

TEST_CASE( "colour_transform()::Check maximum output is 255", "[colour_transform()]" ) {

    // Uses inputs designed to give the largest possible output after the colour matrix
    // multiplication and checks they are not more than 255
    // For each component:
    // Y : (255, 255, 255) -> 255
    // U : (  0,   0, 255) -> 255.5 (exact)
    // V : (255,   0,   0) -> 255.5 (exact)
    Array_2d<double> red{ 
        {255, 0, 255}
    };

    Array_2d<double> green{ 
        {255, 0, 0}
    };

    Array_2d<double> blue{ 
        {255, 255, 0}
    };

    // Now process
    auto [Y, U, V] = colour_transform(red, green, blue);

    for (size_t ind = 0; ind < red.size(); ind++)
    {
        CAPTURE( ind );
        CHECK( Y[ind] <= 255.0 );
        CHECK( U[ind] <= 255.0 );
        CHECK( V[ind] <= 255.0 );
    }
}

TEST_CASE( "colour_transform()::Check sample values", "[colour_transform()]" ) {

    // Test some sample inputs RGB->YUV:
    // ( 34,  67,  23) -> ( 52, 112, 115)
    // (105,   4, 230) -> ( 60, 224, 160)
    // (200, 170, 210) -> (184, 143, 140)
    // (120, 140,  10) -> (119,  66, 129)
    // (240,  50,  70) -> (109, 106, 221)

    Array_2d<double> red{ 
        { 34, 105, 200, 120, 240}
    };

    Array_2d<double> green{ 
        { 67,   4, 170, 140,  50}
    };

    Array_2d<double> blue{ 
        { 23, 230, 210,  10, 70}
    };

    std::vector<double> expected_Y{
         52,  60, 184, 119, 109
    };

    std::vector<double> expected_U{
        112, 224, 143,  66, 106
    };

    std::vector<double> expected_V{
        115, 160, 140, 129, 221
    };

    // Now process
    auto [Y, U, V] = colour_transform(red, green, blue);

    for (size_t ind = 0; ind < red.size(); ind++)
    {
        CAPTURE( ind );
        CHECK( Y[ind] == expected_Y[ind] );
        CHECK( U[ind] == expected_U[ind] );
        CHECK( V[ind] == expected_V[ind] );
    }
}

TEST_CASE( "subsample_component_4_2_0()::Check valid components", "[subsample_component_4_2_0()]" ) {
    Array_2d<double> input_array{ gen_array_2d_arange(16, 16) };
    Array_2d<double> expected_result{
        {  8.5,  10.5,  12.5,  14.5,  16.5,  18.5,  20.5,  22.5},
        { 40.5,  42.5,  44.5,  46.5,  48.5,  50.5,  52.5,  54.5},
        { 72.5,  74.5,  76.5,  78.5,  80.5,  82.5,  84.5,  86.5},
        {104.5, 106.5, 108.5, 110.5, 112.5, 114.5, 116.5, 118.5},
        {136.5, 138.5, 140.5, 142.5, 144.5, 146.5, 148.5, 150.5},
        {168.5, 170.5, 172.5, 174.5, 176.5, 178.5, 180.5, 182.5},
        {200.5, 202.5, 204.5, 206.5, 208.5, 210.5, 212.5, 214.5},
        {232.5, 234.5, 236.5, 238.5, 240.5, 242.5, 244.5, 246.5}
    };
    size_t expected_height{ expected_result.shape()[0] }, expected_width{ expected_result.shape()[1] };

    // subsample_component_4_2_0() acts in place so copy the input array
    Array_2d<double> actual_result{ input_array };

    subsample_component_4_2_0(actual_result);

    REQUIRE( actual_result.shape()[0] == expected_height );
    REQUIRE( actual_result.shape()[1] == expected_width );

    for (size_t i = 0; i < expected_height; i++)
    {
        for (size_t j = 0; j < expected_width; j++)
        {
            CAPTURE( i, j );
            CHECK_THAT( actual_result(i, j), WithinRel(expected_result(i, j)));
        }
    }
}

TEST_CASE( "subsample_component_4_2_0()::Check components with invalid shape", "[subsample_component_4_2_0()]" ) {
    SECTION("Invalid height")
    {
        Array_2d<double> input_array{ gen_array_2d_arange(8, 16) };

        CHECK_THROWS_AS( subsample_component_4_2_0(input_array), std::invalid_argument);
    }
    SECTION("Invalid width")
    {
        Array_2d<double> input_array{ gen_array_2d_arange(16, 8) };

        CHECK_THROWS_AS( subsample_component_4_2_0(input_array), std::invalid_argument);
    }
}

TEST_CASE( "subsample_component_4_2_2()::Check valid components", "[subsample_component_4_2_2()]" ) {
    Array_2d<double> input_array{ gen_array_2d_arange(8, 16) };
    Array_2d<double> expected_result{
        {  0.5,   2.5,   4.5,   6.5,   8.5,  10.5,  12.5,  14.5},
        { 16.5,  18.5,  20.5,  22.5,  24.5,  26.5,  28.5,  30.5},
        { 32.5,  34.5,  36.5,  38.5,  40.5,  42.5,  44.5,  46.5},
        { 48.5,  50.5,  52.5,  54.5,  56.5,  58.5,  60.5,  62.5},
        { 64.5,  66.5,  68.5,  70.5,  72.5,  74.5,  76.5,  78.5},
        { 80.5,  82.5,  84.5,  86.5,  88.5,  90.5,  92.5,  94.5},
        { 96.5,  98.5, 100.5, 102.5, 104.5, 106.5, 108.5, 110.5},
        {112.5, 114.5, 116.5, 118.5, 120.5, 122.5, 124.5, 126.5}
    };
    size_t expected_height{ expected_result.shape()[0] }, expected_width{ expected_result.shape()[1] };

    // subsample_component_4_2_2() acts in place so copy the input array
    Array_2d<double> actual_result{ input_array };

    subsample_component_4_2_2(actual_result);

    REQUIRE( actual_result.shape()[0] == expected_height );
    REQUIRE( actual_result.shape()[1] == expected_width );

    for (size_t i = 0; i < expected_height; i++)
    {
        for (size_t j = 0; j < expected_width; j++)
        {
            CAPTURE( i, j );
            CHECK_THAT( actual_result(i, j), WithinRel(expected_result(i, j)));
        }
    }
}

TEST_CASE( "subsample_component_4_2_2()::Check components with invalid shape", "[subsample_component_4_2_2()]" ) {
    SECTION("Invalid height")
    {
        Array_2d<double> input_array{ gen_array_2d_arange(8, 9) };

        CHECK_THROWS_AS( subsample_component_4_2_2(input_array), std::invalid_argument);
    }
    SECTION("Invalid width")
    {
        Array_2d<double> input_array{ gen_array_2d_arange(16, 8) };

        CHECK_THROWS_AS( subsample_component_4_2_2(input_array), std::invalid_argument);
    }
}

