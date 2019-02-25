#include "bezier.h"
#include <iostream>
#include <chrono>
#include <random>
#include <Eigen/Cholesky>
#include <qpOASES.hpp>
#include "catch.hpp"

TEST_CASE( "example test" ) {
  using namespace splx;
  using namespace std;
  Bezier<double, 2> bez(7.0);

  Bezier<double, 2>::VectorDIM vec;


  int ss = 8;
  while(ss--) {
    vec(0) = vec(1) = 0;
    bez.m_controlPoints.push_back(vec);
  }


  Bezier<double, 2>::QPMatrices QP = bez.getQPMatrices();

  bez.extendQPIntegratedSquaredDerivative(QP, 1, 50);
  bez.extendQPIntegratedSquaredDerivative(QP, 2, 50);


  vec(0) = 7.0;
  vec(1) = -2.0;
  bez.extendQPBeginningConstraint(QP, 0, vec);
  vec(0) = 3;
  vec(1) = 3;
  bez.extendQPBeginningConstraint(QP, 1, vec);
  vec(0) = 0.2;
  vec(1) = 0.5;
  bez.extendQPBeginningConstraint(QP, 2 , vec);

  vec(0) = 3.5;
  vec(1) = 3.5;
  bez.extendQPPositionAt(QP, 3.0, vec, 50.0);
  vec(0) = -3.0;
  vec(1) = 15.0;
  bez.extendQPPositionAt(QP, 4.0, vec, 111100.0);

  vec(0) = -3.5;
  vec(1) = -3.5;
  bez.extendQPPositionAt(QP, 6.0, vec, 50.0);
  vec(0) = -8.5;
  vec(1) = -8.5;
  bez.extendQPPositionAt(QP, 7.0, vec, 10000.0);
  //bez.extendQPDecisionConstraint(QP, -10, 100);



  qpOASES::QProblem qp(QP.x.rows(), QP.A.rows(), qpOASES::HST_SEMIDEF);
  qpOASES::Options options;
  options.setToDefault();
  options.enableRegularisation = qpOASES::BT_FALSE;
  options.enableCholeskyRefactorisation = 0;
  options.printLevel = qpOASES::PL_NONE;
  qp.setOptions(options);
  qpOASES::int_t nWSR = 10000;
  Eigen::LLT<Bezier<double, 2>::Matrix> lltofH(QP.H);




  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> H = QP.H.cast<double>();
  Eigen::Matrix<double, Eigen::Dynamic, 1> g = QP.g.cast<double>();
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A = QP.A.cast<double>();
  Eigen::Matrix<double, Eigen::Dynamic, 1> lbX = QP.lbX.cast<double>();
  Eigen::Matrix<double, Eigen::Dynamic, 1> ubX = QP.ubX.cast<double>();
  Eigen::Matrix<double, Eigen::Dynamic, 1> lbA = QP.lbA.cast<double>();
  Eigen::Matrix<double, Eigen::Dynamic, 1> ubA = QP.ubA.cast<double>();
  Eigen::Matrix<double, Eigen::Dynamic, 1> x = QP.x.cast<double>();
  qpOASES::returnValue status = qp.init(
    H.data(),
    g.data(),
    A.data(),
    lbX.data(),
    ubX.data(),
    lbA.data(),
    ubA.data(),
    nWSR,
    NULL,
    x.data());

    qpOASES::int_t simpleStatus = qpOASES::getSimpleStatus(status);

   REQUIRE( simpleStatus == 0 );




}