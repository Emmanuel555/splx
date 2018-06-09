cmake_minimum_required(VERSION 2.6)
project(bspline)

add_library(bspline bspline.cpp)
add_executable(main main.cpp)

target_link_libraries(main bspline)
