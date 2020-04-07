#ifndef SPLX_PARAMETRIC_CURVE_HPP
#define SPLX_PARAMETRIC_CURVE_HPP
#include <Eigen/Dense>
#include <cmath>
#include <vector>

namespace splx {

template <typename T, unsigned int DIM>
class ParametricCurve {
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


  ParametricCurve(CurveType t) : type(t) {
  }

  virtual ~ParametricCurve() {

  }

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
  * Returns true if the curve is in the non-positive side of the hyperplane
  */
  virtual bool onNonPositiveSide(const Hyperplane& hp) const = 0;

  /*
  * Get i^th control point
  */
  virtual VectorDIM& operator[](std::size_t i) = 0;
  virtual const VectorDIM& operator[](std::size_t i) const = 0;

  /*
  * Get number of control points
  */
  virtual std::size_t numControlPoints() const = 0;


  /*
    Append a new control points to the curve.
  */
  virtual void appendControlPoint(const VectorDIM& cpt) = 0;

  /*
    Remove i^th control point
  */
  virtual void removeControlPoint(std::size_t i) = 0;

  /*
  * Returns the max parameter of the curve, i.e. returns a where curve is defined
  * for u \in [0, a]
  */
  virtual T maxParameter() const = 0;

  /*
  * Sets the new max parameter
  */
  virtual void maxParameter(T nw) = 0;
  
};


}
#endif
