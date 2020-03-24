cmake_minimum_required(VERSION 3.1)
project(splx_tests)

# add_subdirectory(test/Catch2)

file(
GLOB 
all_sources
"test/*.cpp"
)

add_executable(
    splx_tests
    ${all_sources}
)

target_link_libraries(
    splx_tests PUBLIC
    splx
)