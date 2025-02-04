#ifndef JPEG_CPP_JPEG_H
#define JPEG_CPP_JPEG_H

#include <algorithm>
#include <vector>

#include "jpeg_cpp/array.h"
#include "jpeg_cpp/component.h"
#include "jpeg_cpp/encoding.h"
#include "jpeg_cpp/huff_table.h"
#include "jpeg_cpp/general.h"
#include "jpeg_cpp/q_table.h"


namespace JPEG
{

    /// @brief Encodes a greyscale image
    /// @param array_2d Array containing the image data
    /// @param qf Quality factor used for the quantization tables. Defaults to 50, equivalent to spec tables
    /// @return Image data encoded as a JPEG
    std::vector<unsigned char> encode_greyscale_image(const Array_2d<double>& array_2d, int qf=50);

    /// @brief Encodes a colour image
    /// @param red Red component of image
    /// @param green Green component of image
    /// @param blue Blue component of image
    /// @param qf Quality factor used for the quantization tables. Defaults to 50, equivalent to spec tables
    /// @param ss Type of chromiance subsampling to be performed. Defaults to 4:4:4, i.e. no subsampling
    /// @return Image data encoded as a JPEG
    std::vector<unsigned char> encode_colour_image(
        const Array_2d<double>& red, const Array_2d<double>& green, const Array_2d<double>& blue, int qf=50,
        Subsampling ss=Subsampling::ss_4_4_4
    );
}

#endif
