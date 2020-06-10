#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <splx/curve/Bezier.hpp>
#include <splx/internal/bezier.hpp>
#include <stdexcept>
#include <splx/opt/BezierQPOperations.hpp>
#include <splx/opt/PiecewiseCurveQPGenerator.hpp>


TEST_CASE("construction test", "[PiecewiseCurveQPGenerator]") {
    splx::PiecewiseCurveQPGenerator<double, 3> generator;
}