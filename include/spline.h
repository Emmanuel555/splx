#ifndef SPLX_SPLINE_H
#define SPLX_SPLINE_H
#include <Eigen/Dense>

namespace splx {

  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vec;

  class Spline {
  protected:
    unsigned int fac(unsigned int n) const;
    unsigned int comb(unsigned int n, unsigned int k) const;
    unsigned int perm(unsigned int n, unsigned int k) const;
  public:
    /**
     * Evaluate spline at u
    */
    virtual Vec eval(double u) = 0;
    /**
     * Evaluate k^{th} derivative of the spline at u
    */
    virtual Vec eval(double u, unsigned int k) = 0;
  };


}
#endif
