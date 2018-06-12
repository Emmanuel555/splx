#include <qpOASES.hpp>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <random>
#include <iostream>
#include "bspline.h"
#include <chrono>

using std::cout;
using std::endl;

int main() {
  std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_real_distribution<double> distribution(-10, 10);
  auto rand = std::bind(distribution, generator);

  std::vector<splx::Vec> cpts;

  cpts.resize(20);

  for(size_t i = 0; i < cpts.size(); i++) {
    splx::Vec v(2);
    v(0) = i;
    v(1) = (i+rand());
    cpts[i] = v;
    cout << "c" << endl;
    cout << v << endl;
  }

  splx::BSpline bspl(3, 2, 0, 1.0, cpts);

  int i = 0;

  for(double t = 0; t<=1.0; t+=0.001) {
    cout << "d" << endl;
    cout << bspl.eval(t, 0) << endl;
    if((i++)%100 == 0) {
      cout << "v" << endl;
      cout << bspl.eval(t, 1) << endl;
    }
  }
  cout << "d" << endl;
  cout << bspl.eval(1.0, 0) << endl;
  cout << "v" << endl;
  cout << bspl.eval(1.0, 1) << endl;

  return 0;
}
