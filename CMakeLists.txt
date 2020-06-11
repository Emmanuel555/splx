cmake_minimum_required(VERSION 3.9)
project(splx)

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_STANDARD 17)

add_subdirectory(third_party/qp_wrappers)

if(NOT TARGET absl::strings)
    add_subdirectory(third_party/abseil-cpp)
endif()

add_library(splx INTERFACE)
target_link_libraries(
    splx
    INTERFACE
    qp_wrappers_problem
    absl::strings
)
target_include_directories(
    splx INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
 )

include(test/Tests.cmake)