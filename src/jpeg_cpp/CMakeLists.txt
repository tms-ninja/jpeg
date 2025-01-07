
# Add the library
add_library(jpeg_cpp SHARED bit_string.cpp "${jpeg_SOURCE_DIR}/include/jpeg_cpp/bit_string.h")

# Add headers library users need to use this library
target_include_directories(jpeg_cpp PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
                                            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# Enable all warnings including additioal warnings
target_compile_options(jpeg_cpp PRIVATE -Wall -Wextra)

# Ensure compiler supports at least C++17
target_compile_features(jpeg_cpp PUBLIC cxx_std_17)

