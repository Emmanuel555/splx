#ifndef SPLX_SPLINE_H
#define SPLX_SPLINE_H
#include <Eigen/Dense>

namespace splx {

  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vec;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix;

  class Spline {
  protected:
    unsigned int fac(unsigned int n) const;
    unsigned int comb(unsigned int n, unsigned int k) const;
    unsigned int perm(unsigned int n, unsigned int k) const;
  public:
    /**
     * Evaluate k^{th} derivative of the spline at u
    */
    virtual Vec eval(double u, unsigned int k) const = 0;

    /**
     * Return hessian matrix for QP where every element is 0.0.
    */
    virtual Matrix getZeroHessian() const = 0;

    /**
     * Return the g vector for QP where every element is 0.0.
    */
    virtual Vec getZeroG() const = 0;

    /**
     * Converts hessian matrix to upper triangular matrix that can be used in QPs
     * where objective has the term p'Dp instead of 1/2 p'Hp
    */
    virtual Matrix convertHessianToUpperTriangular(const Matrix& H) const;
    virtual void convertHessianToUpperTriangular(Matrix& H) const;
  };


}
#endif
