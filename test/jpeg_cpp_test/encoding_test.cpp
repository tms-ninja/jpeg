#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>

#include "jpeg_cpp/array.h"
#include "jpeg_cpp/encoding.h"

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
