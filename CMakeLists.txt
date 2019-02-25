cmake_minimum_required(VERSION 3.1)
project(splx)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE Debug)

add_compile_options(-std=c++14)
add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/qpOASES)
set(QPOASES_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/third_party/qpOASES/include)

find_package(Eigen3 REQUIRED)

include(CTest)
if (BUILD_TESTING)
      add_subdirectory(test)
endif()
