
# Add the library
add_library(jpeg_cpp SHARED 
            "${jpeg_SOURCE_DIR}/include/jpeg_cpp/array.h"
            array.cpp
            "${jpeg_SOURCE_DIR}/include/jpeg_cpp/bit_string.h"
            bit_string.cpp
            "${jpeg_SOURCE_DIR}/include/jpeg_cpp/component.h"
            component.cpp
            "${jpeg_SOURCE_DIR}/include/jpeg_cpp/encoding.h"
            encoding.cpp
            "${jpeg_SOURCE_DIR}/include/jpeg_cpp/general.h"
            "${jpeg_SOURCE_DIR}/include/jpeg_cpp/huff_table.h"
            huff_table.cpp
            "${jpeg_SOURCE_DIR}/include/jpeg_cpp/jpeg.h"
            jpeg.cpp
            "${jpeg_SOURCE_DIR}/include/jpeg_cpp/q_table.h"
            q_table.cpp
)

# Add headers library users need to use this library
target_include_directories(jpeg_cpp PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
                                            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# Enable all warnings including additioal warnings
if(MSVC)
    target_compile_options(jpeg_cpp PRIVATE /W4)
else()
    target_compile_options(jpeg_cpp PRIVATE -Wall -Wextra)
endif()

# Ensure compiler supports at least C++17
target_compile_features(jpeg_cpp PUBLIC cxx_std_17)

# Only install if installing the python wrapper
if(DEFINED SKBUILD)
    # Think this install the headers, maybe shouldn't go here but in root CMakeLists.txt?
    install(DIRECTORY "${PROJECT_SOURCE_DIR}/include/"
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

    install(
        TARGETS jpeg_cpp
        EXPORT jpeg_cpp_targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

    # # I think the rest of this is something to do with making find_package() work?
    # Think we can skip it for now
    # set(JPEG_CPP_CMAKE_PACKAGE_INSTALL_SUBDIR "share/jpeg_cpp/cmake")

    # # No idea
    # install(
    #     EXPORT jpeg_cpp_targets
    #     NAMESPACE jpeg_cpp::
    #     DESTINATION ${JPEG_CPP_CMAKE_PACKAGE_INSTALL_SUBDIR})

    # include(CMakePackageConfigHelpers)

    # write_basic_package_version_file(
    #     jpeg_cppConfigVersion.cmake
    #     VERSION ${PROJECT_VERSION}
    #     COMPATIBILITY SameMinorVersion)
    
    # configure_package_config_file(
    #     "${PROJECT_SOURCE_DIR}/cmake/jpeg_cppConfig.cmake.in" jpeg_cppConfig.cmake
    #     INSTALL_DESTINATION ${JPEG_CPP_CMAKE_PACKAGE_INSTALL_SUBDIR})
    
    # install(FILES "${PROJECT_BINARY_DIR}/jpeg_cppConfig.cmake"
    #                 "${PROJECT_BINARY_DIR}/jpeg_cppConfigVersion.cmake"
    #         DESTINATION ${JPEG_CPP_CMAKE_PACKAGE_INSTALL_SUBDIR})

endif()
