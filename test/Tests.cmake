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
generate_test(PiecewiseCurveTest)
generate_test(PiecewiseCurveQPGeneratorTest)