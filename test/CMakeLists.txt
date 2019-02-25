cmake_minimum_required(VERSION 3.1)

# if splx is not defined, stop cmaking.
if(NOT TARGET splx)
  message(FATAL_ERROR "splx not defined, do not cmake this file directly")
endif()

# add qpoases without tests
set(QPOASES_BUILD_EXAMPLES CACHE BOOL OFF)
add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/qpOASES qpOASES)
set(QPOASES_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/third_party/qpOASES/include)

# add catch without tests
set(CATCH_BUILD_TESTING CACHE BOOL OFF)
add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/Catch2 Catch2)

# compile main as a seperate library so that we can link it to other
# tests fast
add_library(
  splx_catchmain
  main.cpp
)
target_link_libraries(
  splx_catchmain PUBLIC
  Catch2
)

# these are the libraries that tests link to in general
set(
  TEST_LIBRARIES
  splx_catchmain
  splx
  qpOASES
)

# tests for bezier.h
add_executable(
  splx_bezier_test
  bezier_tests.cpp
)
target_link_libraries(
  splx_bezier_test PRIVATE
  ${TEST_LIBRARIES}
)
target_include_directories(
  splx_bezier_test PRIVATE
  ${QPOASES_INCLUDE_DIR}
)
add_test(NAME splx_bezier_test COMMAND splx_bezier_test)


# tests for bspline.h
add_executable(
  splx_bspline_test
  bspline_tests.cpp
)
target_link_libraries(
  splx_bspline_test PRIVATE
  ${TEST_LIBRARIES}
)
target_include_directories(
  splx_bspline_test PRIVATE
  ${QPOASES_INCLUDE_DIR}
)
add_test(NAME splx_bspline_test COMMAND splx_bspline_test)


# tests for curve.h
add_executable(
  splx_curve_test
  curve_tests.cpp
)
target_link_libraries(
  splx_curve_test PRIVATE
  ${TEST_LIBRARIES}
)
target_include_directories(
  splx_curve_test PRIVATE
  ${QPOASES_INCLUDE_DIR}
)
add_test(NAME splx_curve_test COMMAND splx_curve_test)


# tests for spline.h
add_executable(
  splx_spline_test
  spline_tests.cpp
)
target_link_libraries(
  splx_spline_test PRIVATE
  ${TEST_LIBRARIES}
)
target_include_directories(
  splx_spline_test PRIVATE
  ${QPOASES_INCLUDE_DIR}
)
add_test(NAME splx_spline_test COMMAND splx_spline_test)
