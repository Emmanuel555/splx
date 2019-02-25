cmake_minimum_required(VERSION 3.1)
project(bspline)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE Debug)

add_compile_options(-std=c++14)
add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/qpOASES)
set(QPOASES_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/third_party/qpOASES/include)

find_package(Eigen3 REQUIRED)


include_directories(
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_SOURCE_DIR}/third_party/Catch2/single_include/catch2
  ${EIGEN3_INCLUDE_DIR}
  ${QPOASES_INCLUDE_DIR}

)


add_executable(
  beziertests
  test/bezier_tests.cpp
  test/main.cpp
)

add_executable(
bsplinetests
test/bspline_tests.cpp
test/main.cpp
)


add_executable(
curvetests
test/curve_tests.cpp
test/main.cpp
)


add_executable(
splinetests
test/spline_tests.cpp
test/main.cpp
)



include(CTest)
add_test(NAME BezierTests COMMAND $<TARGET_FILE:beziertests>)
add_test(NAME BSplineTests COMMAND $<TARGET_FILE:bsplinetests>)
add_test(NAME CurveTests COMMAND $<TARGET_FILE:curvetests>)
add_test(NAME SplineTests COMMAND $<TARGET_FILE:splinetests>)
