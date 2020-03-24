#ifndef SPLX_CURVE_H
#define SPLX_CURVE_H
#include <Eigen/Dense>
#include <cmath>
#include <vector>

namespace splx {

template <typename T, unsigned int DIM>
class ParametricCurve {
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

  using VectorDIM = Eigen::Matrix<T, DIM, 1>;
  using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;
  using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using Hyperplane = Eigen::Hyperplane<T, DIM>;
  using MatrixDIM = Eigen::Matrix<T, DIM, DIM, Eigen::RowMajor>;
  using Row = Eigen::Matrix<T, 1, Eigen::Dynamic>;
  using AlignedBox = Eigen::AlignedBox<T, DIM>;
  using ControlPoints = std::vector<VectorDIM, Eigen::aligned_allocator<VectorDIM>>;

  enum class CurveType {
    LINEAR,
    BSPLINE,
    BEZIER
  };

  CurveType type;


  Curve(CurveType t) : type(t) {
  }

  virtual ~Curve() {

  }
  /*
   * The row of k^th derivative of basis functions evaluated at u
   *
   * Remember that basis functions get multiplied with control points to evaluate the function
  */
  virtual Row getBasisRow(T u, unsigned int k) const = 0;

  /**
   * Evaluate k^{th} derivative of the curve at u
  */
  virtual VectorDIM eval(T u, unsigned int k) const = 0;
  
  /*
  * Returns true if the curve is on the negative side of the hyperplane
  * otherwise returns false
  */
  virtual bool onNegativeSide(const Hyperplane& hp) const = 0;

  /*
  * Get i^th control point
  */
  virtual VectorDIM& operator[](std::size_t i) = 0;

  /*
  * Get number of control points
  */
  virtual std::size_t numControlPoints() const = 0;


  /*
  * Returns the max parameter of the curve, i.e. returns a where curve is defined
  * for u \in [0, a]
  */
  virtual T maxParameter() const = 0;
  
};


}
#endif
