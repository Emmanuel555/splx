cmake_minimum_required(VERSION 2.6)
project(bspline)


find_package(Eigen3 REQUIRED)

include_directories(
  ${PROJECT_SOURCE_DIR}/include
  ${EIGEN3_INCLUDE_DIR}
)

add_library(splx src/bspline.cpp)

add_executable(main test/main.cpp)
target_link_libraries(main splx)
