cmake_minimum_required(VERSION 3.1)
project(splx_tests)

function(generate_test testnamebase)
add_executable(
    ${testnamebase}
    test/${testnamebase}.cpp
)

target_link_libraries(
    ${testnamebase} PUBLIC
    splx
)
endfunction(generate_test)

generate_test(BezierTest)
generate_test(BezierQPGeneratorTest)
generate_test(PiecewiseCurveTest)