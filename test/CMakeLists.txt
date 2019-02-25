include_directories(
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_SOURCE_DIR}/third_party/Catch2/single_include/catch2
  ${EIGEN3_INCLUDE_DIR}
  ${QPOASES_INCLUDE_DIR}

)


add_executable(
  beziertests
  bezier_tests.cpp
  main.cpp
)

add_executable(
bsplinetests
bspline_tests.cpp
main.cpp
)


add_executable(
curvetests
curve_tests.cpp
main.cpp
)


add_executable(
splinetests
spline_tests.cpp
main.cpp
)



add_test(NAME BezierTests COMMAND $<TARGET_FILE:beziertests>)
add_test(NAME BSplineTests COMMAND $<TARGET_FILE:bsplinetests>)
add_test(NAME CurveTests COMMAND $<TARGET_FILE:curvetests>)
add_test(NAME SplineTests COMMAND $<TARGET_FILE:splinetests>)