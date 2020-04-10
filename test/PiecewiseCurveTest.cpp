#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include <splx/curve/Bezier.hpp>
#include <splx/curve/PiecewiseCurve.hpp>
#include <Eigen/Dense>
#include <iostream>

TEST_CASE("", "[PiecewiseCurve]") {

    using VectorDIM = Eigen::Matrix<double, 3, 1>;

    double double_eq_epsilon = 1e-13;

    splx::PiecewiseCurve<double, 3> piecewiseCurve;

    {
        splx::Bezier<double, 3> piece1 {3.5, {{1, 2, 3}, {3, 2, 1}, {2, 2.3, 3.4}, {3, 2.2, 3.1}, {2.2314, 2.231, 
    2.22}, {-2.11, -.231, 1.2}, {-5, -5, -5}}};
        splx::Bezier<double, 3> piece2 {2, {{-5, -5, -5}, {3, 2, 1}, {2, 2.3, 3.4}, {3, 2.2, 3.1}, {2.2314, 
    2.231, 2.22}, {-2.11, -.231, 1.2}, {-5, -5, -5}}};
    
        piecewiseCurve.addPiece(piece1);
        piecewiseCurve.addPiece(piece2);
    }


    REQUIRE(piecewiseCurve.numPieces() == 2);
    REQUIRE(piecewiseCurve.type(0) == splx::ParametricCurve<double, 3>::CurveType::BEZIER);
    REQUIRE(piecewiseCurve.type(1) == splx::ParametricCurve<double, 3>::CurveType::BEZIER);
    REQUIRE_THROWS_AS(piecewiseCurve.type(2), std::domain_error);

    REQUIRE((piecewiseCurve.eval(3.4, 0) - piecewiseCurve[0].eval(3.4, 0)).squaredNorm() < double_eq_epsilon);
    REQUIRE((piecewiseCurve.eval(3.7, 0) - piecewiseCurve[1].eval(0.2, 0)).squaredNorm() < double_eq_epsilon);

    REQUIRE((piecewiseCurve.eval(3.4, 1) - piecewiseCurve[0].eval(3.4, 1)).squaredNorm() < double_eq_epsilon);
    REQUIRE((piecewiseCurve.eval(3.7, 1) - piecewiseCurve[1].eval(0.2, 1)).squaredNorm() < double_eq_epsilon);

    REQUIRE_THROWS_AS(piecewiseCurve.eval(5.50000001, 2), std::domain_error);


    REQUIRE((piecewiseCurve.eval(2.123616263123, 5) - VectorDIM{9.831787628214174, 2.5970117330035167, -9.036818965112952}).squaredNorm() < double_eq_epsilon);
    REQUIRE((piecewiseCurve.eval(4.91273162413, 6) - VectorDIM{-133.52625000000006, 37.698749999999976, -10.12499999999995}).squaredNorm() < double_eq_epsilon);
}