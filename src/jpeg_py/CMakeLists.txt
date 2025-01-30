# Adapted from the example at https://github.com/scikit-build/scikit-build-sample-projects/tree/main/projects/hello-cmake-package

find_package(
  Python
  COMPONENTS Interpreter Development.Module NumPy
  REQUIRED
)

python_add_library(jpeg MODULE jpeg.cpp WITH_SOABI)

# Link the the C++ header of jpeg_cpp
target_link_libraries(jpeg PRIVATE jpeg_cpp PUBLIC Python::NumPy)

if(MSVC)
    target_compile_options(jpeg_cpp PRIVATE /W4)
else()
    target_compile_options(jpeg_cpp PRIVATE -Wall -Wextra)
endif()

# Only install if initialized from scikit-build-core
if(DEFINED SKBUILD)
    install(TARGETS jpeg DESTINATION .)

    # Ensure the C++ library is loaded when the Python extension is loaded
    if(WIN32)
        # On Windows, copy the resulting DLL to the same directory as the Python wrapper (.pyd)
        add_custom_command(
            TARGET jpeg POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_RUNTIME_DLLS:jpeg> ${SKBUILD_PLATLIB_DIR}
            COMMAND_EXPAND_LISTS
        )
    else()
        # On Unix, can use rpath to tell it where to find the share library
        if(APPLE)
            set_target_properties(
                jpeg PROPERTIES INSTALL_RPATH "@loader_path/${CMAKE_INSTALL_LIBDIR}")
        else()
            set_target_properties(jpeg PROPERTIES INSTALL_RPATH
                                                    "$ORIGIN/${CMAKE_INSTALL_LIBDIR}")
        endif()
    endif()
endif()
