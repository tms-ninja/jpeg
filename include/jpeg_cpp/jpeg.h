#ifndef JPEG_CPP_JPEG_H
#define JPEG_CPP_JPEG_H

#include <vector>

#include "jpeg_cpp/array.h"
#include "jpeg_cpp/component.h"
#include "jpeg_cpp/encoding.h"
#include "jpeg_cpp/huff_table.h"
#include "jpeg_cpp/general.h"
#include "jpeg_cpp/q_table.h"


namespace JPEG
{
    std::vector<unsigned char> encode_greyscale_image(const Array_2d<double>& array_2d);
}

#endif
