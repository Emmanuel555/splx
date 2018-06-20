cmake_minimum_required(VERSION 3.1)
project(bspline)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Debug)

find_package(Eigen3 REQUIRED)

include_directories(
  ${PROJECT_SOURCE_DIR}/include
  ${EIGEN3_INCLUDE_DIR}
)

add_library(
splx
src/spline.cpp
src/bspline.cpp
)

add_executable(splxtest test/main.cpp)
target_link_libraries(splxtest splx)
