# Adapted from the example at https://github.com/scikit-build/scikit-build-sample-projects/tree/main/projects/hello-cmake-package

find_package(
  Python
  COMPONENTS Interpreter Development.Module NumPy
  REQUIRED
)

python_add_library(jpeg MODULE jpeg.cpp WITH_SOABI)

# Link the the C++ header of jpeg_cpp
target_link_libraries(jpeg PRIVATE jpeg_cpp PUBLIC Python::NumPy)

# Only install if initialized from scikit-build-core
if(DEFINED SKBUILD)
    install(TARGETS jpeg DESTINATION .)

    # Something to do with ensuring the C++ library is loaded when the Python extension is loaded
    if(APPLE)
        set_target_properties(
            jpeg PROPERTIES INSTALL_RPATH "@loader_path/${CMAKE_INSTALL_LIBDIR}")
    else()
        set_target_properties(jpeg PROPERTIES INSTALL_RPATH
                                                  "$ORIGIN/${CMAKE_INSTALL_LIBDIR}")
    endif()
endif()
