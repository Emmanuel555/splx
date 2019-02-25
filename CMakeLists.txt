cmake_minimum_required(VERSION 3.1)
project(bspline)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -pedantic ")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE Debug)

add_compile_options(-std=c++14)
add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/qpOASES)
set(QPOASES_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/third_party/qpOASES/include)

find_package(Eigen3 REQUIRED)


include_directories(
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_SOURCE_DIR}/third_party
  ${EIGEN3_INCLUDE_DIR}
  ${QPOASES_INCLUDE_DIR}

)


add_executable(
  catchbezier
  test/catch_bezier.cpp
  test/main.cpp
)

target_link_libraries(
catchbezier
qpOASES
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


include(CTest)
add_test(NAME RunTests COMMAND $<TARGET_FILE:catchbezier>)