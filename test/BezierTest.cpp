#include "catch.hpp"
#include <splx/Bezier.hpp>
#include <stdexcept>

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

        SECTION("getBasisRow") {
            for(unsigned int k = 0; k < 10; k++) {
                Row row = bez.getBasisRow(0, k);
                REQUIRE(row.cols() == 0);
            }

            REQUIRE_THROWS_AS(bez.getBasisRow(-1, 0), std::domain_error);
            REQUIRE_THROWS_AS(bez.getBasisRow(0.0001, 0), std::domain_error);
        }

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