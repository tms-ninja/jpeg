cmake_minimum_required(VERSION 3.28)

project(
  jpeg
  VERSION 0.1
  LANGUAGES CXX)

# By default, no symbols are exported on Windows unless they are marked __declspec(dllexport)
# This ensures all symbols are exported and accessible from the Python wrapper
set( CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS True )

include(GNUInstallDirs)

# Add directory containing the source code for the factorial_cpp library
add_subdirectory(src)

# Add test
include(CTest)

# Only build tests if we're not building the python module 
if(NOT DEFINED SKBUILD)
    add_subdirectory(test)
endif()
