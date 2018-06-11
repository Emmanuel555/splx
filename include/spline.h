#ifndef SPLX_SPLINE_H
#define SPLX_SPLINE_H
#include <Eigen/Dense>

namespace splx {

  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vec;

  class Spline {
  public:
    /**
     * Evaluate spline at u
    */
    virtual Vec eval(double u) = 0;
    /**
     * Evaluate n^{th} derivative of the spline at u
    */
    virtual Vec eval(double u, int n) = 0;
  };


}
#endif
