cmake_minimum_required(VERSION 3.28)

project(
  jpeg
  VERSION 0.1
  LANGUAGES CXX)

include(GNUInstallDirs)


# Add directory containing the source code for the factorial_cpp library
add_subdirectory(src)

# Add test
include(CTest)

# Only build tests if we're not building the python module 
if(NOT DEFINED SKBUILD)
    add_subdirectory(test)
endif()
