#include <qpOASES.hpp>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <random>
#include <iostream>
#include "bspline.h"
#include <chrono>
#include <qpOASES.hpp>
#include <Eigen/Cholesky>

using std::cout;
using std::endl;

int main() {
  std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_real_distribution<double> distribution(-10, 10);
  auto rand = std::bind(distribution, generator);

  std::vector<splx::Vec> cpts;


  splx::Vec vec(2);
  vec(0) = 0.0;
  vec(1) = 0.0;
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);
  cpts.push_back(vec);

  splx::BSpline bspl(4, 2, 0, 7.0, cpts);
  bspl.generateUniformKnotVector();


  splx::QPMatrices QP = bspl.getQPMatrices();

  //bspl.printKnotVector();

  bspl.extendQPIntegratedSquaredDerivative(QP, 1, 1.0);
  bspl.extendQPIntegratedSquaredDerivative(QP, 2, 1.0);
  vec(0) = -7.0;
  vec(1) = -7.0;
  bspl.extendQPBeginningConstraint(QP, 0, vec);
  vec(0) = 3;
  vec(1) = 3;
  bspl.extendQPBeginningConstraint(QP, 1, vec);
  vec(0) = 0.2;
  vec(1) = 0.5;
//  bspl.extendQPBeginningConstraint(QP, 2 , vec);

  vec(0) = 3.5;
  vec(1) = 3.5;
  bspl.extendQPPositionAt(QP, 3.0, vec, 50.0);
  vec(0) = 10.0;
  vec(1) = 2.0;
  bspl.extendQPPositionAt(QP, 4.0, vec, 100.0);

  vec(0) = -3.5;
  vec(1) = -3.5;
  bspl.extendQPPositionAt(QP, 6.0, vec, 50.0);
  vec(0) = -8.5;
  vec(1) = -8.5;
  bspl.extendQPPositionAt(QP, 7.0, vec, 1000.0);
  bspl.extendQPDecisionConstraint(QP, -10, 10);

  splx::Hyperplane hp(2);
  hp.normal()(0) = 2;
  hp.normal()(1) = 1;
  hp.offset() = 0;

  bspl.extendQPHyperplaneConstraint(QP, 5U, 25U, hp);


  qpOASES::QProblem qp(QP.x.rows(), QP.A.rows(), qpOASES::HST_SEMIDEF);
  qpOASES::Options options;
  options.setToMPC();
  options.printLevel = qpOASES::PL_NONE;
  qp.setOptions(options);
  qpOASES::int_t nWSR = 10000;

  std::cerr << "H" << std::endl << QP.H << std::endl << std::endl;
  std::cerr << "g" << std::endl << QP.g << std::endl << std::endl;
  std::cerr << "A" << std::endl << QP.A << std::endl << std::endl;
  std::cerr << "lbX" << std::endl << QP.lbX << std::endl << std::endl;
  std::cerr << "ubX" << std::endl << QP.ubX << std::endl << std::endl;
  std::cerr << "lbA" << std::endl << QP.lbA << std::endl << std::endl;
  std::cerr << "ubA" << std::endl << QP.ubA << std::endl << std::endl;
  std::cerr << "x" << std::endl << QP.x << std::endl << std::endl;

  qpOASES::returnValue status = qp.init(
    QP.H.data(),
    QP.g.data(),
    QP.A.data(),
    QP.lbX.data(),
    QP.ubX.data(),
    QP.lbA.data(),
    QP.ubA.data(),
    nWSR,
    NULL,
    QP.x.data());

  qpOASES::int_t simpleStatus = qpOASES::getSimpleStatus(status);

  if(simpleStatus != 0) {
    std::cerr << std::endl << "failed " << std::endl;
  }

  std::cerr << "status: " << simpleStatus << std::endl;

  qp.getPrimalSolution(QP.x.data());

  std::cerr << "x optimized" << std::endl << QP.x << std::endl;

  bspl.loadControlPoints(QP);

  for(int i=0; i < bspl.m_controlPoints.size(); i++) {
    cout << "c" << endl;
    cout << bspl.m_controlPoints[i] << endl;
  }

  int i = 0;

  for(double t = bspl.m_a; t<bspl.m_b; t+=0.01) {
    cout << "d" << endl;
    cout << bspl.eval(t, 0) << endl;
    if((i++) % 50 == 0) {
      cout << "v" << endl;
      cout << bspl.eval(t, 1) << endl;
    }
  }
  cout << "d" << endl;
  cout << bspl.eval(bspl.m_b, 0) << endl;
  cout << "v" << endl;
  cout << bspl.eval(bspl.m_b, 1) << endl;
  return 0;
}
