cmake_minimum_required(VERSION 3.1)
project(splx)

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_STANDARD 17)

add_subdirectory(third_party/qp_wrappers)
add_subdirectory(third_party/abseil-cpp)

find_package(Eigen3 REQUIRED)

add_library(splx INTERFACE)
target_link_libraries(
splx 
INTERFACE
qp_wrappers
absl::strings
)
target_include_directories(
  splx INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${EIGEN3_INCLUDE_DIR}
 )

include(test/Tests.cmake)