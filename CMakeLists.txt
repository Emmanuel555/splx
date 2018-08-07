cmake_minimum_required(VERSION 3.1)
project(bspline)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Debug)

find_package(Eigen3 REQUIRED)


include_directories(
  ${PROJECT_SOURCE_DIR}/include
  ${EIGEN3_INCLUDE_DIR}
  ${QPOASES_INCLUDE_DIR}
)


add_executable(
bsplinetest
test/testbspline.cpp
)

target_link_libraries(
bsplinetest
qpOASES
)

add_executable(
beziertest
test/testbezier.cpp
)

target_link_libraries(
beziertest
qpOASES
)

add_executable(
splinetest
test/testspline.cpp
)

target_link_libraries(
splinetest
qpOASES
)
