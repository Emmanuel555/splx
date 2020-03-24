#ifndef SPLX_BEZIER_H
#define SPLX_BEZIER_H

#include "curve.h"
#include <iostream>
#include <numeric>
#include <Eigen/StdVector>
#include <algorithm>
#include <memory>
#include <stdexcept>


namespace splx {

template<typename T, unsigned int DIM>
class Bezier : public Curve<T, DIM> {
  public:
    using VectorDIM = typename Curve<T, DIM>::VectorDIM;
    using Vector = typename Curve<T, DIM>::Vector;
    using Hyperplane = typename Curve<T, DIM>::Hyperplane;
    using Matrix = typename Curve<T, DIM>::Matrix;
    using Row = typename Curve<T, DIM>::Row;
    using ControlPoints = typename Curve<T, DIM>::ControlPoints;
    using CurveType = typename Curve<T, DIM>::CurveType;

    Bezier(T a, const ControlPoints& cpts) : Curve<T, DIM>(Curve<T, DIM>::CurveType::BEZIER),
                                                      m_controlPoints(cpts),
                                                      m_a(a)
    {
      if(m_a <= 0) {
        throw std::domain_error(
          "max parameter should be positive. given "
          + std::to_string(m_a)
        );
      }
    }

    Bezier(T a) : Curve<T, DIM>(Curve<T, DIM>::CurveType::BEZIER), m_a(a) {
      if(m_a <= 0) {
        throw std::domain_error(
          "max parameter should be positive. given "
          + std::to_string(m_a)
        );
      }
    }

    /*
    * Construct from curve pointer
    * curve->m_type should be beizer
    */
    Bezier(const std::shared_ptr<Curve<T, DIM>> curve): Curve<T, DIM>(Curve<T, DIM>::CurveType::BEZIER, curve->a) {
      if(curve->type != CurveType::Bezier) {
        throw std::domain_error(
          "Tried to initialize bezier curve with another type of curve"s
        );
      }

      auto bezpt = std::static_pointer_cast<Bezier<T, DIM>>(curve);
      m_controlPoints = bezpt->m_controlPoints;
      m_a = bezpt->m_a;
    }

    ~Bezier() {

    }

    Bezier(const Bezier<T, DIM>& rhs) = delete;
    Bezier(Bezier<T, DIM>&& rhs) = delete;

    Bezier& operator=(const Bezier<T, DIM>& rhs) = delete;
    Bezier& operator=(Bezier<T, DIM>&& rhs) = delete;

    std::size_t numControlPoints() const override {
      return m_controlPoints.size();
    }

    VectorDIM& operator[](std::size_t i) override {
      return m_controlPoints[i];
    }

    T maxParameter() const override {
      return m_a;
    }

    unsigned int degree() const {
      return this->numControlPoints - 1;
    }

    Row getQPBasisRow(T u, unsigned int k) const override {
      if(k < 0) {
        throw std::domain_error(
          "can't take k^th derivative of the curve since k is negative. k = "s
          + std::to_string(k)
        );
      }

      if(u < 0 || u > maxParameter()) {
        throw std::domain_error(
          "u is outside of the range [0, "s
          + std::to_string(maxParameter)
          + "]"s
        );
      }

      unsigned int degree = this->degree();
      Row result(degree + 1);
      T oneOverA = 1/m_a;
      for(unsigned int i = 0; i <= degree; i++) {
        T base = 0.0;
        T mult = 1.0;
        for(unsigned int j = 0; j+k <= degree; j++, mult *= u) {
          if(j+k >= i) {
            base += pow(oneOverA, i) * this->comb(degree-i, j+k-i)
                  * pow(-oneOverA, j+k-i) * this->perm(j+k, k)
                  * mult;
          }
        }
        base *= this->comb(degree, i);
        result(i) = base;
      }
      return result;
    }


    /*
     * Evaluate k^th derivative of bezier curve at u = u.
     *
     * @fails if u > m_a
    */
    VectorDIM eval(T u, unsigned int k) const override {
      if(k < 0) {
        throw std::domain_error(
          "can't take k^th derivative of the curve since k is negative. k = "s
          + std::to_string(k)
        );
      }

      if(u < 0 || u > maxParameter()) {
        throw std::domain_error(
          "u is outside of the range [0, "s
          + std::to_string(maxParameter)
          + "]"s
        );
      }

      Row basis = this->getQPBasisRow(u, k);
      VectorDIM result;
      result.setZero();

      for(typename std::size_t i = 0; i < this->numControlPoints(); i++) {
        result += (*this)[i] * basis(i);
      }

      return result;
    }


    /*
    * Returns true if the curve is in the negative side of the hyperplane hp
    * leverage convex hull property
    */
    bool onNegativeSide(const Hyperplane& hp) const override {
      for(const auto& cp: m_controlPoints) {
        if(hp.signedDistance(cp) > 0)
          return false;
      }
      return true;
    }

  private:
    /**
     * Bezier curve is defined for u \in [0, m_a]
    */
    T m_a;

    /**
     * Control points of the bezier curve.
     * Implicitly defines the order of the curve
     *
     * General equation of bezier curves is
     * f(u | m_a) = \sum_{i=0}^{d} P_i C(d, i) (u/m_a)^i (1-u/m_a)^(d-i)
     * Therefore, the degree d of the bezier is m_controlPoints.size() - 1
    */
    ControlPoints m_controlPoints;
};

}

#endif
