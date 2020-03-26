#ifndef SPLX_BEZIER_HPP
#define SPLX_BEZIER_HPP

#include <splx/ParametricCurve.hpp>
#include <iostream>
#include <numeric>
#include <Eigen/StdVector>
#include <algorithm>
#include <memory>
#include <stdexcept>


namespace splx {

template<typename T, unsigned int DIM>
class Bezier : public ParametricCurve<T, DIM> {
  public:
    using Base = ParametricCurve<T, DIM>;

    using VectorDIM = typename ParametricCurve<T, DIM>::VectorDIM;
    using Vector = typename ParametricCurve<T, DIM>::Vector;
    using Hyperplane = typename ParametricCurve<T, DIM>::Hyperplane;
    using Matrix = typename ParametricCurve<T, DIM>::Matrix;
    using Row = typename ParametricCurve<T, DIM>::Row;
    using ControlPoints = typename ParametricCurve<T, DIM>::ControlPoints;
    using CurveType = typename ParametricCurve<T, DIM>::CurveType;

    Bezier(): Base(CurveType::BEZIER), m_a(0) {

    }

    Bezier(T a, const ControlPoints& cpts) : Base(CurveType::BEZIER),
                                                      m_a(a),
                                                      m_controlPoints(cpts)
                                                      
    {
      if(m_a < 0) {
        throw std::domain_error(
          std::string("max parameter should be non-negative. given ")
          + std::to_string(m_a)
        );
      }
    }

    Bezier(T a) : Base(CurveType::BEZIER), m_a(a) {
      if(m_a < 0) {
        throw std::domain_error(
          std::string("max parameter should be non-negative. given ")
          + std::to_string(m_a)
        );
      }
    }

    /*
    * Construct from curve pointer
    * curve->m_type should be beizer
    */
    Bezier(const std::shared_ptr<ParametricCurve<T, DIM>> curve): Base(CurveType::BEZIER, curve->a) {
      if(curve->type != CurveType::Bezier) {
        throw std::domain_error(
          std::string("tried to initialize bezier curve with another type of curve")
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
      return this->m_controlPoints[i];
    }

    const VectorDIM& operator[](std::size_t i) const override {
      return this->m_controlPoints[i];
    }

    void appendControlPoint(const VectorDIM& cpt) override {
      this->m_controlPoints.push_back(cpt);
    }

    void removeControlPoint(std::size_t idx) override {
      if(idx >= this->m_controlPoints.size()) {
        throw std::domain_error(
          std::string("given index out of range for removal")
        );
      }

      this->m_controlPoints.erase(this->m_controlPoints.begin() + idx);
    }

    T maxParameter() const override {
      return m_a;
    }

    void maxParameter(T nw) override {
      if(nw < 0) {
        throw std::domain_error(
          std::string("max parameter should be non-negative. given ")
          + std::to_string(nw)
        );
      }

      m_a = nw;
    }

    unsigned int degree() const {
      return this->numControlPoints() - 1;
    }

    Row getBasisRow(T u, unsigned int k) const override {
      if(u < 0 || u > maxParameter()) {
        throw std::domain_error(
          std::string("u is outside of the range [0, ")
          + std::to_string(maxParameter())
          + std::string("]")
        );
      }

      if(m_a == 0) {
        Row result(this->numControlPoints());
        result.setZero();
        if(k == 0 && this->numControlPoints() > 0) {
          result(0) = 1.0;
        }
        return result;
      }

      unsigned int degree = this->degree();
      Row result(degree + 1);

      T oneOverA = 1/m_a;
      for(unsigned int i = 0; i <= degree; i++) {
        T base = 0.0;
        T mult = 1.0;
        for(unsigned int j = 0; j+k <= degree; j++, mult *= u) {
          if(j+k >= i) {
            // base += pow(oneOverA, i) * this->comb(degree-i, j+k-i)
            //       * pow(-oneOverA, j+k-i) * this->perm(j+k, k)
            //       * mult;

            base += this->comb(degree-i, j+k-i)
              * pow(oneOverA, j+k) * this->perm(j+k, k)
              * mult * ((j+k-i)%2 == 0 ? 1 : -1);
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
      if(u < 0 || u > maxParameter()) {
        throw std::domain_error(
          std::string("u is outside of the range [0, ")
          + std::to_string(maxParameter())
          + std::string("]")
        );
      }

      Row basis = this->getBasisRow(u, k);
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
        if(hp.signedDistance(cp) >= 0)
          return false;
      }

      if(this->numControlPoints() == 0) {
        VectorDIM zeros;
        zeros.setZero();
        return hp.signedDistance(zeros) < 0;
      }

      return true;
    }


    /*
    * Returns true if the curve is in the negative side of the hyperplane hp
    * leverage convex hull property
    */
    bool onNonPositiveSide(const Hyperplane& hp) const override {
      for(const auto& cp: m_controlPoints) {
        if(hp.signedDistance(cp) > 0)
          return false;
      }

      if(this->numControlPoints() == 0) {
        VectorDIM zeros;
        zeros.setZero();
        return hp.signedDistance(zeros) <= 0;
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
