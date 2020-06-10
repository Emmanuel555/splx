#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <splx/curve/Bezier.hpp>
#include <splx/internal/bezier.hpp>
#include <stdexcept>
#include <splx/opt/BezierQPOperations.hpp>


TEST_CASE("initialization of Bezier curves", "[bezier]") {
    using Bez = splx::Bezier<double, 3>;
    using VectorDIM = Bez::VectorDIM;
    using Row = Bez::Row;
    using Hyperplane = Bez::Hyperplane;
    using ControlPoints = Bez::ControlPoints;

    VectorDIM zero_vecdim;
    zero_vecdim.setZero();

    double double_eq_epsilon = 1e-13;

    SECTION("bezier curve with no control points") {
        splx::Bezier<double, 3> bez;


        REQUIRE(bez.maxParameter() == 0);
        REQUIRE(bez.numControlPoints() == 0);

        SECTION("eval") {
            for(unsigned int k = 0; k < 10; k++) {
                REQUIRE(bez.eval(0, 0) == zero_vecdim);
            }

            REQUIRE_THROWS_AS(bez.eval(0.0001, 0), std::domain_error);
            REQUIRE_THROWS_AS(bez.eval(-1, 0), std::domain_error);
        }

        SECTION("onNegativeSide") {
            VectorDIM normal(1, 1, 1);
            Hyperplane hp(normal, -1);

            REQUIRE(bez.onNegativeSide(hp));

            Hyperplane hp2(normal, 0);
            REQUIRE(!bez.onNegativeSide(hp2));
        }

        SECTION("onNonPositiveSide") {
            VectorDIM normal(1, 1, 1);
            Hyperplane hp(normal, -1);
            REQUIRE(bez.onNonPositiveSide(hp));

            Hyperplane hp2(normal, 0);
            REQUIRE(bez.onNonPositiveSide(hp2));
        }
    }

    SECTION("bezier curve *with* a lot of control points") {
        SECTION("linear bezier curve") {
            ControlPoints cpts;
            for(int i = 0; i < 100; i++) {
                VectorDIM cpt(i, i, i);
                cpts.push_back(cpt);
            }

            Bez bez(111, cpts);

            REQUIRE(bez.maxParameter() == 111);
            bez.maxParameter(120);
            REQUIRE(bez.maxParameter() == 120);
            REQUIRE(bez.numControlPoints() == 100);

            for(int i = 0; i < 100; i++) {
                REQUIRE(bez[i] == VectorDIM(i, i, i));
            }

            REQUIRE(bez.degree() == 99);

            SECTION("onNegativeSide & onPositiveSide") {
                VectorDIM normal(1,1,1);
                Hyperplane hp(normal, 4);
                REQUIRE(!bez.onNegativeSide(hp));
                REQUIRE(!bez.onNonPositiveSide(hp));

                VectorDIM normal2(-1, -1, 2);
                Hyperplane hp2(normal2, -1);
                REQUIRE(bez.onNegativeSide(hp2));
                REQUIRE(bez.onNonPositiveSide(hp2));


                Hyperplane hp3(normal2, 0);
                REQUIRE(!bez.onNegativeSide(hp3));
                REQUIRE(bez.onNonPositiveSide(hp3));
            }
        }
    }

    SECTION("bezier with control points") {
        Bez bez;

        bez.appendControlPoint(VectorDIM(1, 2, 3));
        bez.appendControlPoint(VectorDIM(2, 2, 3));
        bez.appendControlPoint(VectorDIM(3, 2, 1));

        REQUIRE(bez.numControlPoints() == 3);

        bez.removeControlPoint(1);
        REQUIRE(bez[0] == VectorDIM(1, 2, 3));
        REQUIRE(bez[1] == VectorDIM(3, 2, 1));
        REQUIRE(bez.numControlPoints() == 2);

        bez.appendControlPoint(VectorDIM(2, 2.3, 3.4));
        bez.appendControlPoint(VectorDIM(3, 2.2, 3.1));
        bez.appendControlPoint(VectorDIM(2.2314, 2.231, 2.22));
        bez.appendControlPoint(VectorDIM(-2.11, -.231, 1.2));
        bez.appendControlPoint(VectorDIM(-5, -5, -5));

        REQUIRE(bez.numControlPoints() == 7);

        bez.maxParameter(3.5);
        REQUIRE(bez.maxParameter() == 3.5);

        REQUIRE_THROWS_AS(bez.eval(-1, 0), std::domain_error);
        REQUIRE_THROWS_AS(bez.eval(3.50000001, 0), std::domain_error);

        REQUIRE(bez.eval(1.2, 9) == zero_vecdim);
    
        REQUIRE((bez.eval(0.1, 0) - VectorDIM(1.3083683839412488, 2.0033587482309585, 2.7077993750767115)).norm() < double_eq_epsilon);
        REQUIRE((bez.eval(0, 0) - bez[0]).norm() < double_eq_epsilon);
        REQUIRE((bez.eval(3.5, 0) - bez[6]).norm() < double_eq_epsilon);
        REQUIRE((bez.eval(0.111 ,0) - VectorDIM(1.3383694241843003, 2.0040973541448315, 2.6814193887593163)).norm() < double_eq_epsilon);
        REQUIRE((bez.eval(0.22, 1) - VectorDIM(2.1149983786062885, 0.11889184332111918, -1.5010478716032567)).norm() < double_eq_epsilon);
        REQUIRE((bez.eval(0.37, 2) - VectorDIM(-3.486583859387949, 0.16434710152038, 4.813030336097375)).norm() < double_eq_epsilon);
        REQUIRE((bez.eval(0.12, 3) - VectorDIM(11.575138623483749, -1.6469262980642423, -17.32425004901019)).norm() < double_eq_epsilon);
        REQUIRE((bez.eval(2.76, 6) - VectorDIM(-2.2987319909221484, 4.054212785489039, 2.7808821154450962)).norm() < double_eq_epsilon);
        REQUIRE((bez.eval(2.76, 6) - bez.eval(3.4, 6)).norm() < double_eq_epsilon);
        REQUIRE((bez.eval(2.22, 7) - zero_vecdim).norm() < double_eq_epsilon);

    }
}


TEST_CASE("splx::internal::bezier::bernsteinCoefficientMatrix test") {
    double double_eq_epsilon = 1e-15;

    splx::BezierQPGenerator<double, 3> generator(8, 0.35914595895879409);

    auto bern = splx::internal::bezier::bernsteinCoefficientMatrix(7, 0.35914595895879409, 1);
    Eigen::Matrix<double, 8, 8> mtr;
    mtr << -19.490682897543422,19.490682897543422,0,0,0,0,0,0,
            325.61718841079289,-651.2343768215859,325.61718841079289,0,0,0,0,0,
            -2266.6076304658627,6799.822891397589,-6799.8228913975881,2266.6076304658627,0,0,0,0,
            8414.8058244880795,-33659.223297952318,50488.834946928473,-33659.223297952318,8414.8058244880795,0,0,0,
            -17572.53342530342,87862.667126517103,-175725.33425303421,175725.33425303421,-87862.667126517103,17572.53342530342,0,0,
            19571.467239946949,-117428.80343968168,293572.00859920419,-391429.34479893895,293572.00859920419,-117428.80343968165,19571.467239946949,0,
            -9082.4109955586973,63576.876968910881,-190730.63090673267,317884.38484455441,-317884.38484455441,190730.63090673267,-63576.876968910881,9082.4109955586973,
            0,0,0,0,0,0,0,0;
    REQUIRE((bern.transpose()-mtr).squaredNorm() < double_eq_epsilon);

    bern = splx::internal::bezier::bernsteinCoefficientMatrix(7, 0.35914595895879409, 2);
    mtr << 325.61718841079289,-651.2343768215859,325.61718841079289,0,0,0,0,0,-4533.2152609317254,13599.645782795178,-13599.645782795176,4533.2152609317254,0,0,0,0,25244.417473464237,-100977.66989385695,151466.5048407854,-100977.66989385695,25244.417473464237,0,0,0,-70290.133701213679,351450.66850606841,-702901.33701213682,702901.33701213682,-351450.66850606841,70290.133701213679,0,0,97857.336199734738,-587144.01719840837,1467860.0429960208,-1957146.7239946947,1467860.0429960208,-587144.01719840826,97857.336199734738,0,-54494.465973352191,381461.26181346527,-1144383.7854403961,1907306.3090673264,-1907306.3090673264,1144383.7854403961,-381461.26181346527,54494.465973352191,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    REQUIRE((bern.transpose()-mtr).squaredNorm() < double_eq_epsilon);

    bern = splx::internal::bezier::bernsteinCoefficientMatrix(7, 0.35914595895879409, 4);
    mtr << 50488.834946928473,-201955.33978771389,302933.00968157081,-201955.33978771389,50488.834946928473,0,0,0,-421740.8022072821,2108704.0110364105,-4217408.0220728209,4217408.0220728209,-2108704.0110364105,421740.8022072821,0,0,1174288.0343968167,-7045728.2063809,17614320.515952252,-23485760.687936336,17614320.515952252,-7045728.2063808991,1174288.0343968167,0,-1089889.3194670437,7629225.2362693055,-22887675.708807919,38146126.181346528,-38146126.181346528,22887675.708807919,-7629225.2362693055,1089889.3194670437,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;

    REQUIRE((bern.transpose()-mtr).squaredNorm() < double_eq_epsilon);
}
