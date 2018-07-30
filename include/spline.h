#ifndef SPLX_SPLINE_H
#define SPLX_SPLINE_H
#include <Eigen/Dense>

namespace splx {

template <typename T, unsigned int DIM>
class Spline {
protected:
  unsigned int fac(unsigned int n) const {
    unsigned int res = 1;
    for(unsigned int i = 2; i<=n; i++)
      res *= i;
    return res;
  }
  unsigned int comb(unsigned int n, unsigned int k) const {
    k = std::min(k, n-k);
    unsigned int top = 1;
    unsigned int bottom = 1;
    for(unsigned int i=0; i<k; i++) {
      bottom *= (i+1);
      top *= (n-i);
    }
    return top / bottom;
  }
  unsigned int perm(unsigned int n, unsigned int k) const {
    return fac(n) / fac(n-k);
  }
public:

  using VectorDIM = Eigen::Matrix<T, DIM, 1>;
  using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;
  using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using Hyperplane = Eigen::Hyperplane<T, DIM>;

  /**
  * QP is assumed to be formulated as 1/2x^THx + x^Tg
  * subject to
  * lb <= Ax <= ub
  * where x is a column vector of decision variables.
  */
  typedef struct QPMatrices {
    Matrix H;
    Vector g;
    Matrix A;
    Vector lbA;
    Vector ubA;
    Vector lbX;
    Vector ubX;
    Vector x;
  } QPMatrices;

  /**
   * Evaluate k^{th} derivative of the spline at u
  */
  virtual VectorDIM eval(double u, unsigned int k) const = 0;


  /**
   * Every element of H is 0.0
   * Every element of g is 0.0
   * Every element of x is the current values of decision variables.
  */
  virtual QPMatrices getQPMatrices() const = 0;


  /**
   * Converts hessian matrix to upper triangular matrix that can be used in QPs
   * where objective has the term p'Dp instead of 1/2 p'Hp
  */
  virtual Matrix convertHessianToUpperTriangular(const Matrix& H) const {
    unsigned int S = H.rows();
    Matrix U(H);
    for(unsigned int i = 0; i < S; i++) {
      U(i, i) /= 2.0;
      for(unsigned int j = 0; j < i; j++) {
        U(i, j) = 0.0;
      }
    }
    return U;
  }
  virtual void convertHessianToUpperTriangular(Matrix& H) const {
    unsigned int S = H.rows();
    for(unsigned int i = 0; i < S; i++) {
      H(i, i) /= 2.0;
      for(unsigned int j = 0; j < i; j++) {
        H(i, j) = 0.0;
      }
    }
  }
};


}
#endif
