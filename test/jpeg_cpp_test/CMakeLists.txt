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

if(MSVC)
    target_compile_options(jpeg_cpp PRIVATE /W4)
else()
    target_compile_options(jpeg_cpp PRIVATE -Wall -Wextra)
endif()

# Copy the jpeg_cpp DLL to the same directory as the test executable, otherwise Windows won't
# be able to find it
if(WIN32)
    add_custom_command(
        TARGET jpeg_cpp_test POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_RUNTIME_DLLS:jpeg_cpp_test> $<TARGET_FILE_DIR:jpeg_cpp_test>
        COMMAND_EXPAND_LISTS
    )
endif()

# Add the actual test
# Need this for Catch to find the tests
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)

include(Catch)
catch_discover_tests(jpeg_cpp_test)
