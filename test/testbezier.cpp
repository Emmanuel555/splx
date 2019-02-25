#include "bezier.h"
#include <iostream>
#include <chrono>
#include <random>
#include <Eigen/Cholesky>
#include <qpOASES.hpp>

int main() {
  using namespace splx;
  using namespace std;


//  std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
//  std::uniform_real_distribution<double> distribution(-10, 10);
//  auto rand = std::bind(distribution, generator);

  Bezier<double, 2> bez(7.0);

  Bezier<double, 2>::VectorDIM vec;



  /*vec(0) = -1.9331079999999932;
  vec(1) = -6.999999999999952e-06;
  bez.m_controlPoints.push_back(vec);
  vec(0) = -1.9052911428571366;
  vec(1) = -9.8571428571428e-06;
  bez.m_controlPoints.push_back(vec);
  vec(0) = -1.870046047619042;
  vec(1) = -1.3428571428571362e-05;
  bez.m_controlPoints.push_back(vec);
  vec(0) = -1.8274856285714225;
  vec(1) = -1.7657142857142786e-05;
  bez.m_controlPoints.push_back(vec);
  vec(0) = -1.77836785714285;
  vec(1) = -2.2399999999999914e-05;
  bez.m_controlPoints.push_back(vec);
  vec(0) = -1.723203190476181;
  vec(1) = -2.7523809523809403e-05;
  bez.m_controlPoints.push_back(vec);
  vec(0) = -1.6624898571428424;
  vec(1) = -3.2857142857142674e-05;
  bez.m_controlPoints.push_back(vec);
  vec(0) = -1.5966739999999755;
  vec(1) = -3.7999999999999704e-05;
  bez.m_controlPoints.push_back(vec);*/

  /*int ss = 10;
  while(ss--) {
    vec(0) = rand();
    vec(1) = rand();
    bez.m_controlPoints.push_back(vec);
  }*/

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

  std::cerr << "H" << std::endl << QP.H << std::endl << std::endl;
  std::cerr << "g" << std::endl << QP.g << std::endl << std::endl;
  std::cerr << "A" << std::endl << QP.A << std::endl << std::endl;
  std::cerr << "lbX" << std::endl << QP.lbX << std::endl << std::endl;
  std::cerr << "ubX" << std::endl << QP.ubX << std::endl << std::endl;
  std::cerr << "lbA" << std::endl << QP.lbA << std::endl << std::endl;
  std::cerr << "ubA" << std::endl << QP.ubA << std::endl << std::endl;
  std::cerr << "x" << std::endl << QP.x << std::endl << std::endl;
  Eigen::LLT<Bezier<double, 2>::Matrix> lltofH(QP.H);
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

  bez.loadControlPoints(QP);

  for(unsigned int i = 0; i < bez.m_controlPoints.size(); i++) {
    cout << "c" << endl;
    cout << bez.m_controlPoints[i] << endl;
  }


  int i = 0;
  for(double t = 0; t < bez.m_a; t += 0.01, i++) {
    cout << "d" << endl << bez.eval(t, 0) << endl;
    if(i%51 == 0) {
      cout << "v" << endl << bez.eval(t, 1) << endl;
    }
  }
  cout << "d" << endl << bez.eval(bez.m_a, 0) << endl;
  cout << "v" << endl << bez.eval(bez.m_a, 1) << endl;

  return 0;
}
