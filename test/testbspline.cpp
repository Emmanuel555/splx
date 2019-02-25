#include <qpOASES.hpp>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <random>
#include <iostream>
#include "bspline.h"
#include <chrono>
#include <qpOASES.hpp>
#include <Eigen/Cholesky>
#include <Eigen/StdVector>

using std::cout;
using std::endl;

int main() {
  using namespace splx;

  //std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
  //std::uniform_real_distribution<double> distribution(-10, 10);
  //auto rand = std::bind(distribution, generator);

  std::vector<BSpline<double, 2>::VectorDIM, Eigen::aligned_allocator<BSpline<double, 2>::VectorDIM> > cpts;


  BSpline<double, 2>::VectorDIM vec;
  vec(0) = 0.0;
  vec(1) = 0.0;
  int i = 20;
  while(i--) {
    cpts.push_back(vec);
  }
  BSpline<double, 2> bspl(4, 7.0, cpts);
  bspl.generateClampedUniformKnotVector();
  //bspl.printKnotVectorNumbered();


  BSpline<double, 2>::QPMatrices QP = bspl.getQPMatrices();

  //bspl.printKnotVector();

  bspl.extendQPIntegratedSquaredDerivative(QP, 1, 100.0);
  //bspl.extendQPIntegratedSquaredDerivative(QP, 2, 15.5);
  vec(0) = -2.0;
  vec(1) = -2.0;
  bspl.extendQPBeginningConstraint(QP, 0, vec);
  vec(0) = 3;
  vec(1) = 3;
  bspl.extendQPBeginningConstraint(QP, 1, vec);
  vec(0) = 0.2;
  vec(1) = 0.5;
  bspl.extendQPBeginningConstraint(QP, 2 , vec);

  vec(0) = 3.5;
  vec(1) = 3.5;
  bspl.extendQPPositionAt(QP, 3.0, vec, 50.0);
  vec(0) = -3.0;
  vec(1) = 15.0;
  bspl.extendQPPositionAt(QP, 4.0, vec, 100.0);

  vec(0) = -3.5;
  vec(1) = -3.5;
  bspl.extendQPPositionAt(QP, 6.0, vec, 50.0);
  vec(0) = 10;
  vec(1) = 20;
  bspl.extendQPPositionAt(QP, 7.0, vec, 10000.0);
  //bspl.extendQPDecisionConstraint(QP, -100, 100);

  BSpline<double, 2>::Hyperplane hp(2);
  hp.normal()(0) = 1;
  hp.normal()(1) = 1;
  hp.offset() = 0;

  bspl.extendQPHyperplaneConstraint(QP, hp);


  qpOASES::QProblem qp(QP.x.rows(), QP.A.rows(), qpOASES::HST_SEMIDEF);
  qpOASES::Options options;
  options.setToDefault();
  options.enableRegularisation = qpOASES::BT_FALSE;
  options.enableCholeskyRefactorisation = 0;
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
  Eigen::LLT<BSpline<double, 2>::Matrix> lltofH(QP.H);
  if(lltofH.info() == Eigen::NumericalIssue) {
    std::cerr << "non psd" << endl;
  } else {
    std::cerr << "psd" << std::endl;
  }
  int a; std::cin >> a;



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

  if(simpleStatus != 0) {
    std::cerr << std::endl << "failed " << std::endl;
    int a;
    std::cin >> a;
  }

  //std::cerr << "status: " << simpleStatus << std::endl;

  qp.getPrimalSolution(x.data());

  QP.x = x;

  std::cerr << "x optimized" << std::endl << QP.x << std::endl;

  bspl.loadControlPoints(QP);

  for(unsigned int i=0; i < bspl.m_controlPoints.size(); i++) {
    cout << "c" << endl;
    cout << bspl.m_controlPoints[i] << endl;
  }


  for(double t = 0; t<bspl.m_a; t+=0.01) {
    cout << "d" << endl;
    cout << bspl.eval(t, 0) << endl;
    if((i++) % 50 == 0) {
      cout << "v" << endl;
      cout << bspl.eval(t, 1) << endl;
    }
  }
  cout << "d" << endl;
  cout << bspl.eval(bspl.m_a, 0) << endl;
  cout << "v" << endl;
  cout << bspl.eval(bspl.m_a, 1) << endl;
  return 0;
}
