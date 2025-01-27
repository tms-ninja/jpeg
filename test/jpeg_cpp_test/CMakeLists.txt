Include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.4.0 # or a later release
)

FetchContent_MakeAvailable(Catch2)

add_executable(jpeg_cpp_test 
              array_test.cpp
              bit_string_test.cpp
              encoding_test.cpp
              huff_table_test.cpp
              jpeg_test.cpp
              q_table_test.cpp
)


target_link_libraries(jpeg_cpp_test PRIVATE jpeg_cpp Catch2::Catch2WithMain)

# Add the actual test
# Need this for Catch to find the tests
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)

include(Catch)
catch_discover_tests(jpeg_cpp_test)
