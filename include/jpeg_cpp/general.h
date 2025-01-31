#ifndef GENERAL_H
#define GENERAL_H

namespace JPEG
{
    /// @brief Represents the type of an image component
    enum class Image_Component
    {
        Luminance,
        Chrominance
    };

    // Width and height of a data unit
    constexpr unsigned int du_width{ 8 };
    constexpr unsigned int du_height{ 8 };
}

#endif
