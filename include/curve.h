#ifndef SPLX_CURVE_H
#define SPLX_CURVE_H
#include <Eigen/Dense>
#include <cmath>

namespace splx {

template <typename T, unsigned int DIM>
class Curve {
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
  T pow(T base, unsigned int exp) const {
    if(base == 0 && exp == 0) {
      return 1;
    }
    return std::pow(base, exp);
  }
public:

  enum CurveType {
    BSPLINE,
    BEZIER
  };

  CurveType m_type;

  Curve(CurveType t) : m_type(t) {

  }

  using VectorDIM = Eigen::Matrix<T, DIM, 1>;
  using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;
  using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using Hyperplane = Eigen::Hyperplane<T, DIM>;
  using MatrixDIM = Eigen::Matrix<T, DIM, DIM, Eigen::RowMajor>;
  using Row = Eigen::Matrix<T, 1, Eigen::Dynamic>;
  using AlignedBox = Eigen::AlignedBox<T, DIM>;

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
   * Evaluate k^{th} derivative of the curve at u
  */
  virtual VectorDIM eval(T u, unsigned int k) const = 0;


  /**
   * Every element of H is 0.0
   * Every element of g is 0.0
   * Every element of x is the current values of decision variables.
  */
  virtual QPMatrices getQPMatrices() const = 0;


  /*
   * Get parameter span
  */
  virtual T parameterSpan() const = 0;

  virtual void extendQPIntegratedSquaredDerivative(QPMatrices& QP, unsigned int k, T lambda) const = 0;
  virtual void extendQPPositionAt(QPMatrices& QP, T u, const VectorDIM& pos, T theta) const = 0;
  virtual void extendQPBeginningConstraint(QPMatrices& QP, unsigned int k, const VectorDIM& target) const = 0;
  virtual void extendQPDecisionConstraint(QPMatrices&QP, T lb, T ub) const = 0;
  virtual void loadControlPoints(const QPMatrices& QP) = 0;
  virtual void extendQPHyperplaneConstraint(QPMatrices& QP, const Hyperplane& hp) const = 0;
  /*
   * The row of k^th derivative of basis functions evaluated at u
   *
   * Remember that basis functions get multiplied with control points to evaluate the function
  */
  virtual Row getQPBasisRow(T u, unsigned int k) const = 0;
  // resets the A, lbA, ubA matrices
  virtual void resetConstraintMatrices(QPMatrices& QP) const = 0;

  // resets H matrix
  virtual void resetCostMatrix(QPMatrices& QP) const = 0;

  /*
  * Returns true if the curve is on the negative side of the hyperplane
  * otherwise returns false
  */
  virtual bool onNegativeSide(const Hyperplane& hp) const = 0;

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

  /*
  * Returns true if the curve intersects with the aligned box at any point from
  * parameter from to parameter to when robot is a box with side 2*radius
  */
  virtual bool intersects(const AlignedBox& box, T from, T to, T radius) const {
    const T step = 0.01;
    const VectorDIM rad(radius, radius, radius);
    for(T t = from; t <= to; t += step) {
      VectorDIM pos = eval(t, 0);
      AlignedBox robot_box(pos - rad, pos + rad);
      if(box.intersects(robot_box)) {
        return true;
      }
    }
    return false;
  }
};


}
#endif
