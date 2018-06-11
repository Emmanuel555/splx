#include <qpOASES.hpp>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <random>
#include <iostream>
#include "bspline.h"

int main() {
  std::default_random_engine generator;
  std::uniform_real_distribution<double> dist(-0.5, 0.5);
  std::vector<splx::Vec> cpts;

  cpts.resize(20);

  for(int i = 0; i < cpts.size(); i++) {
    splx::Vec v(2);
    v(0) = i;
    v(1) = (i+dist(generator)) * 0.2;
    cpts[i] = v;
  }


  splx::BSpline bspl(3, 2, 0, 1.0, cpts);

  bspl.printKnotVectorNumbered();


  return 0;
}
