#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <splx/curve/Bezier.hpp>
#include <splx/internal/bezier.hpp>
#include <stdexcept>
#include <splx/opt/BezierQPOperations.hpp>
#include <splx/opt/PiecewiseCurveQPGenerator.hpp>

using namespace std;

TEST_CASE("construction test", "[PiecewiseCurveQPGenerator]") {
    splx::PiecewiseCurveQPGenerator<double, 3> generator;
    splx::Bezier<double, 3> bez;
    using VectorDIM = Eigen::Matrix<double, 3, 1>;
    vector<int> num_bezier_control_points = {8,8,8,8};

 
    bez.appendControlPoint(VectorDIM(1, 2, 3));
    bez.appendControlPoint(VectorDIM(2, 2, 3));
    bez.appendControlPoint(VectorDIM(3, 2, 1));

    for(const auto& cpts: num_bezier_control_points) {
        cout << cpts << endl; // 8
        generator.addBezier(cpts, 0);
    }
    cout << generator.numPieces() << endl; // 4 
    
}