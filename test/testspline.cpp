#include "spline.h"
#include "bezier.h"
#include <Eigen/Dense>
#include <qpOASES.hpp>
#include <random>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>

std::vector<std::string> split(std::string inp, char delim = ',') {
  using namespace std;
  string cur;
  vector<string> res;
  for(unsigned int i = 0; i < inp.size(); i++) {
    if(inp[i] == delim) {
      res.push_back(cur);
      cur.clear();
    } else {
      cur.push_back(inp[i]);
    }
  }
  res.push_back(cur);
  return res;
}

int main(int argc, char** argv) {
  using namespace std;
  using namespace splx;

/*  if(argc < 2) {
    std::cout << "no arguments" << std::endl;
    return 0;
  }


  string line;

  std::ifstream inp(argv[1]);
  getline(inp, line);


  Spline<double, 2> spl;


  while(std::getline(inp, line)) {
    if(line.size() == 0) {
      continue;
    }


    vector<string> tokens = split(line);

    Bezier<double, 2> bez(stod(tokens[0]));
    for(unsigned int i = 1; i< 16; i+=2) {
      Bezier<double, 2>::VectorDIM vec;
      vec(0) = stod(tokens[i]);
      vec(1) = stod(tokens[i+1]);
    //  cout << vec << endl;
      bez.m_controlPoints.push_back(vec);
    }

    spl.addPiece(bez);
  }

  //int a; cin >> a;

  int i = 0;
  for(double t = 0; t < spl.totalSpan(); t+=0.01) {
    cout << "d" << endl;
    cout << spl.eval(t, 0) << endl;
    if((i++) % 50 == 0) {
      cout << "v" << endl;
      cout << spl.eval(t, 1) << endl;
    }

  }
  cout << "d" << endl;
  cout << spl.eval(spl.totalSpan(), 0) << endl;
  cout << "v" << endl;
  cout << spl.eval(spl.totalSpan(), 1) << endl;


*/

  Spline<double, 2> spl;
  Bezier<double, 2> bez(0.17);
  Bezier<double, 2>::VectorDIM vec;
  vec(0) = 1;
  vec(1) = 2;
  bez.m_controlPoints.push_back(vec);
  vec(0) = 1.2;
  vec(1) = -2.3;
  bez.m_controlPoints.push_back(vec);
  vec(0) = -2.2;
  vec(1) = 1.1;
  bez.m_controlPoints.push_back(vec);
  vec(0) = -2.2;
  vec(1) = 1.1;
  bez.m_controlPoints.push_back(vec);
  vec(0) = -2.2;
  vec(1) = 1.1;
  bez.m_controlPoints.push_back(vec);
  vec(0) = -2.2;
  vec(1) = 1.1;
  bez.m_controlPoints.push_back(vec);
  vec(0) = -2.2;
  vec(1) = 1.1;
  bez.m_controlPoints.push_back(vec);
  vec(0) = -2.2;
  vec(1) = 1.1;
  bez.m_controlPoints.push_back(vec);

  spl.addPiece(bez);

  Bezier<double, 2> bez2(0.32);
  vec(0) = 3;
  vec(1) = 3;
  bez2.m_controlPoints.push_back(vec);
  vec(0) = 2.7;
  vec(1) = 12;
  bez2.m_controlPoints.push_back(vec);
  vec(0) = 2.7;
  vec(1) = 12;
  bez2.m_controlPoints.push_back(vec);
  vec(0) = 2.7;
  vec(1) = 12;
  bez2.m_controlPoints.push_back(vec);
  vec(0) = 2.7;
  vec(1) = 12;
  bez2.m_controlPoints.push_back(vec);
  vec(0) = 2.7;
  vec(1) = 12;
  bez2.m_controlPoints.push_back(vec);

  spl.addPiece(bez2);
  spl.addPiece(bez2);
  spl.addPiece(bez2);
  spl.addPiece(bez2);
  spl.addPiece(bez2);
//  spl.addPiece(bez2);
  //spl.addPiece(bez);

  std::vector<Spline<double, 2>::QPMatrices> qps = spl.getQPMatrices();

  spl.extendQPIntegratedSquaredDerivative(qps, 1, 5.02);
  spl.extendQPIntegratedSquaredDerivative(qps, 2, 100.02);

  vec(0) = 1;
  vec(1) = 0.2;
  spl.extendQPBeginningConstraint(qps, 0, vec);

  vec(0) = -1;
  vec(1) = -2;
  spl.extendQPBeginningConstraint(qps, 1, vec);



  vec(0) = 0;
  vec(1) = 0;
  spl.extendQPPositionAt(qps, spl.totalSpan()-1e-10, vec, 1000.7);

  vec(0) = -3;
  vec(1) = 2;
//  spl.extendQPPositionAt(qps, 8.1, vec, 50.7);

  //spl.extendQPDecisionConstraint(qps, -100, 100);

  Spline<double, 2>::Hyperplane hp;
  hp.normal()(0) = 1;
  hp.normal()(1) = 1;
  hp.offset() = 0;

//  spl.extendQPHyperplaneConstraint(qps, 3, hp);


  Spline<double, 2>::QPMatrices QP = spl.combineQPMatrices(qps);

  spl.extendQPContinuityConstraint(QP, qps, 0, 0);
  spl.extendQPContinuityConstraint(QP, qps, 0, 1);
  spl.extendQPContinuityConstraint(QP, qps, 0, 2);
  spl.extendQPContinuityConstraint(QP, qps, 0, 3);
  spl.extendQPContinuityConstraint(QP, qps, 0, 4);
  spl.extendQPContinuityConstraint(QP, qps, 0, 5);

  spl.extendQPContinuityConstraint(QP, qps, 1, 0);
  spl.extendQPContinuityConstraint(QP, qps, 1, 1);
  spl.extendQPContinuityConstraint(QP, qps, 1, 2);
  spl.extendQPContinuityConstraint(QP, qps, 1, 3);
  spl.extendQPContinuityConstraint(QP, qps, 1, 4);
  spl.extendQPContinuityConstraint(QP, qps, 1, 5);

  spl.extendQPContinuityConstraint(QP, qps, 2, 0);
  spl.extendQPContinuityConstraint(QP, qps, 2, 1);
  spl.extendQPContinuityConstraint(QP, qps, 2, 2);
  spl.extendQPContinuityConstraint(QP, qps, 2, 3);
  spl.extendQPContinuityConstraint(QP, qps, 2, 4);
  spl.extendQPContinuityConstraint(QP, qps, 2, 5);

  spl.extendQPContinuityConstraint(QP, qps, 3, 0);
  spl.extendQPContinuityConstraint(QP, qps, 3, 1);
  spl.extendQPContinuityConstraint(QP, qps, 3, 2);
  spl.extendQPContinuityConstraint(QP, qps, 3, 3);
  spl.extendQPContinuityConstraint(QP, qps, 3, 4);
  spl.extendQPContinuityConstraint(QP, qps, 3, 5);

  spl.extendQPContinuityConstraint(QP, qps, 4, 0);
  spl.extendQPContinuityConstraint(QP, qps, 4, 1);
  spl.extendQPContinuityConstraint(QP, qps, 4, 2);
  spl.extendQPContinuityConstraint(QP, qps, 4, 3);
  spl.extendQPContinuityConstraint(QP, qps, 4, 4);
  spl.extendQPContinuityConstraint(QP, qps, 4, 5);

  //spl.extendQPContinuityConstraint(QP, qps, 0, 1);
  //spl.extendQPContinuityConstraint(QP, qps, 1, 1);


  int a;

  qpOASES::QProblem qp(QP.x.rows(), QP.A.rows(), qpOASES::HST_SEMIDEF);
  qpOASES::Options options;
  options.setToDefault();
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
  std::cin >> a;



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
    std::cin >> a;
  }

  //std::cerr << "status: " << simpleStatus << std::endl;

  qp.getPrimalSolution(x.data());

  QP.x = x;

  std::cerr << "x optimized" << std::endl << QP.x << std::endl;

  std::cerr << "constraint mult" << std::endl << QP.A * QP.x << std::endl;
  cin >> a;
  spl.loadControlPoints(QP, qps);

  std::cerr << spl.totalSpan() << std::endl;

  int i = 0;
  for(double t = 0; t < spl.totalSpan(); t += 0.01, i++) {
    cout << "d" << endl << spl.eval(t, 0) << endl;
    if(i%51 == 0) {
      cout << "v" << endl << spl.eval(t, 1) << endl;
    }
  }
  cout << "d" << endl << spl.eval(spl.totalSpan() - 1e-10, 0) << endl;
  cout << "v" << endl << spl.eval(spl.totalSpan() - 1e-10, 1) << endl;
  return 0;
}
