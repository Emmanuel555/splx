#ifndef SPLX_BEZIER_H
#define SPLX_BEZIER_H

#include "curve.h"
#include <iostream>
#include <numeric>
#include <Eigen/StdVector>
#include <algorithm>

using std::max;

namespace splx {

template<typename T, unsigned int DIM>
class Bezier : public Curve<T, DIM> {
  public:
    using VectorDIM = typename Curve<T, DIM>::VectorDIM;
    using Vector = typename Curve<T, DIM>::Vector;
    using Hyperplane = typename Curve<T, DIM>::Hyperplane;
    using Matrix = typename Curve<T, DIM>::Matrix;
    using QPMatrices = typename Curve<T, DIM>::QPMatrices;
    using Row = typename Curve<T, DIM>::Row;


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
    std::vector<VectorDIM, Eigen::aligned_allocator<VectorDIM> > m_controlPoints;

    /*
     * @fails if a < 0
    */
    Bezier(T a, const std::vector<VectorDIM, Eigen::aligned_allocator<VectorDIM> >& cpts) : Curve<T, DIM>(Curve<T, DIM>::CurveType::BEZIER),
                                                      m_a(a), m_controlPoints(cpts)
                                          {
      assert(a >= 0);
    }

    /*
     * @fails if a < 0
    */
    Bezier(T a) : Curve<T, DIM>(Curve<T, DIM>::CurveType::BEZIER), m_a(a) {
      assert(a >= 0);
    }

    /*
    * Construct from curve pointer
    * curve->m_type should be beizer
    */
    Bezier(const std::shared_ptr<Curve<T, DIM>> curve): Curve<T, DIM>(Curve<T, DIM>::CurveType::BEZIER) {
      assert((curve->m_type == Curve<T, DIM>::CurveType::BEZIER));
      auto bezpt = std::static_pointer_cast<Bezier<T, DIM>>(curve);
      m_a = bezpt->m_a;
      m_controlPoints = bezpt->m_controlPoints;
    }

    /**
     * In all functions that extends these matrices, it is assumed that
     * order of variables is
     * p0x, p1x, ..., pdx, p0y, p1y, ...
    */
    QPMatrices getQPMatrices() const override {
      unsigned int S = DIM * m_controlPoints.size(); // number of decision variables
      QPMatrices QP;
      QP.H.resize(S, S);
      QP.g.resize(S);
      QP.A.resize(0, S);
      QP.lbA.resize(0);
      QP.ubA.resize(0);
      QP.lbX.resize(S);
      QP.ubX.resize(S);
      QP.x.resize(S);

      for(unsigned int i = 0; i < S; i++) {
        for(unsigned int j = 0; j < S; j++) {
          QP.H(i, j) = 0.0;
        }
        QP.g(i) = 0.0;
        QP.lbX(i) = std::numeric_limits<T>::lowest();
        QP.ubX(i) = std::numeric_limits<T>::max();
      }

      for(unsigned int i = 0; i < DIM; i++) {
        for(typename std::vector<VectorDIM, Eigen::aligned_allocator<VectorDIM> >::size_type j = 0; j < m_controlPoints.size(); j++) {
          QP.x(i*m_controlPoints.size() + j) = m_controlPoints[j](i);
        }
      }

      return QP;
    }

    /*
     * Evaluate k^th derivative of bezier curve at u = u.
     *
     * @fails if u > m_a
    */
    VectorDIM eval(T u, unsigned int k) const override {
      assert(u <= m_a);
      std::vector<T> basis = evalBasisFuncs(u, k);
      VectorDIM result;
      for(unsigned int i = 0; i < DIM; i++) {
        result(i) = 0;
      }

      for(typename std::vector<T>::size_type i = 0; i < basis.size(); i++) {
        result += m_controlPoints[i] * basis[i];
      }

      return result;
    }

    /*
     * Add integral from 0 to m_a of the square of norm of k^th derivative of bezier
     * to QP hessian with scale factor lambda
     *
     * order of variables are p0x, p1x, ... pnx, p0y, ..., pny, ...
    */
    void extendQPIntegratedSquaredDerivative(QPMatrices& QP, unsigned int k, T lambda) const override {
      unsigned int degree = m_controlPoints.size() - 1;
      if(k > degree) {
        return;
      }

      Matrix D(degree+1, degree+1); // get k^th derivative coefficients from a0 + a1u + ... a(degree)u^degree
      for(unsigned int m = 0; m < degree+1; m++) {
        for(unsigned int n = 0; n < degree+1; n++) {
          D(m, n) = 0.0;
        }
        if(m <= degree-k) {
          D(m, m+k) = Curve<T, DIM>::perm(m+k, k);
        }
      }

      Matrix M = getBasisCoefficientMatrix().transpose();

      Matrix SQI(degree+1, degree+1);

      for(unsigned int i = 0; i <= degree; i++) {
        for(unsigned int j = 0; j <= degree; j++) {
          SQI(i, j) = 2.0 * pow(m_a, i+j+1) / (i+j+1);
        }
      }

      Matrix Hext = lambda * M.transpose() * D.transpose() * SQI * D * M;
      for(unsigned int d = 0; d < DIM; d++) {
        QP.H.block(d*m_controlPoints.size(), d*m_controlPoints.size(), m_controlPoints.size(), m_controlPoints.size()) += Hext;
      }
    }

    /**
     * Add the cost theta * ||f(u) - pos||^2 to H and g.
     *
     * order of variables is assumed to be
     * p0x, p1x, p2x, ..., pnx, p0y, p1y, ..., pny, ...
    */

    void extendQPPositionAt(QPMatrices& QP, T u, const VectorDIM& pos, T theta) const override {
      assert(u <= m_a);

      std::vector<T> basis = evalBasisFuncs(u, 0);
      Vector M(m_controlPoints.size());
      for(unsigned int i = 0; i<m_controlPoints.size(); i++) {
        M(i) = basis[i];
      }

      Matrix Hext = 2 * theta * M * M.transpose();

      for(unsigned int d = 0; d < DIM; d++) {
        Vector Gext = -2 * theta * pos[d] * M;
        QP.g.block(d*m_controlPoints.size(), 0, m_controlPoints.size(), 1) += Gext;
        QP.H.block(d*m_controlPoints.size(), d*m_controlPoints.size(), m_controlPoints.size(), m_controlPoints.size()) += Hext;
      }
    }

    /**
     * Add constraint that requires the k^th derivative of bezier at u=0 to be target.
    */
    void extendQPBeginningConstraint(QPMatrices& QP, unsigned int k, const VectorDIM& target) const override {
      unsigned int degree = m_controlPoints.size();
      assert(k <= degree);
      unsigned int ridx = QP.A.rows();
      QP.A.conservativeResize(QP.A.rows() + DIM, QP.A.cols());
      QP.lbA.conservativeResize(QP.lbA.rows() + DIM);
      QP.ubA.conservativeResize(QP.ubA.rows() + DIM);

      std::vector<T> basis = evalBasisFuncs(0, k);
      unsigned int S = m_controlPoints.size() * DIM;


      for(unsigned int d = 0; d < DIM; d++) {
        for(unsigned int m = 0; m < S; m++) {
          if(m >= d*m_controlPoints.size() && m < ((d+1)*m_controlPoints.size())) {
            QP.A(ridx+d, m) = basis[m - d*m_controlPoints.size()];
          } else {
            QP.A(ridx+d, m) = 0.0;
          }
        }
        QP.lbA(ridx+d) = target(d);
        QP.ubA(ridx+d) = target(d);
      }
    }

    /**
     * Add constraint that requires curve to be on the negative side of hp
     *
     * Effectively, enforces all points to be in the negative side of the hp.
    */
    void extendQPHyperplaneConstraint(QPMatrices& QP, const Hyperplane& hp) const override {
      unsigned int ridx = QP.A.rows();
      QP.A.conservativeResize(QP.A.rows() + m_controlPoints.size(), QP.A.cols());
      QP.lbA.conservativeResize(QP.lbA.rows() + m_controlPoints.size());
      QP.ubA.conservativeResize(QP.ubA.rows() + m_controlPoints.size());


      unsigned int S = DIM * m_controlPoints.size();

      for(unsigned int i = 0; i < m_controlPoints.size(); i++) {
        for(unsigned int s = 0; s < S; s++) {
          QP.A(ridx + i, s) = 0.0;
        }
        for(unsigned int d = 0; d < DIM; d++) {
          QP.A(ridx + i, d*m_controlPoints.size() + i) = hp.normal()(d);
        }
        QP.lbA(ridx + i) = std::numeric_limits<T>::lowest();
        QP.ubA(ridx + i) = -hp.offset();
      }
    }

    /**
     * Add constraints that requires all decision variables to be less than ub and more than lb.
     *
     * Effectively, it constraints the curve to be inside a box.
    */
    void extendQPDecisionConstraint(QPMatrices& QP, T lb, T ub) const override {
      for(unsigned int i=0; i < QP.lbX.rows(); i++) {
        QP.lbX(i) = lb;
        QP.ubX(i) = ub;
      }
    }

    Row getQPBasisRow(T u, unsigned int k) const override {
      std::vector<T> basis = evalBasisFuncs(u, k);
      assert(basis.size() == m_controlPoints.size());
      Row M(m_controlPoints.size());
      for(unsigned int i = 0; i<m_controlPoints.size(); i++) {
        M(i) = basis[i];
      }

      return M;
    }

    /**
     * Load control points from QP.x
    */
    void loadControlPoints(const QPMatrices& QP) override {
      for(unsigned int i = 0; i < m_controlPoints.size(); i++) {
        for(unsigned int d = 0; d < DIM; d++) {
          m_controlPoints[i](d) = QP.x(d * m_controlPoints.size() + i);
        }
      }
    }

/*
    void resetCostMatrix(QPMatrices& QP) const override {
      unsigned int S = DIM * m_controlPoints.size();
      for(int i = 0; i < S; i++)
        for(int j = 0; j < S; j++)
          QP.H(i, j) = 0.0;
    }

    void resetConstraintMatrices(QPMatrices& QP) const override {
      unsigned int S = DIM * m_controlPoints.size();
      QP.A.resize(0, S);
      QP.lbA.resize(0);
      QP.ubA.resize(0);
    }
*/

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
     * Evaluate k^{th} derivative of basis functions at u = u, where the bezier curve is parametrized by
     * m_a, and degree of basis functions is m_controlPoints.size() - 1
    */
    std::vector<T> evalBasisFuncs(T u, unsigned int k) const {
      unsigned int degree = m_controlPoints.size() - 1;
      std::vector<T> basis(degree + 1, 0.0);
      T oneOverA = 1 / m_a;
      for(unsigned int i = 0; i <= degree; i++) {
        T base = 0.0;
        T mult = 1.0;
        for(unsigned int j = 0; j <= degree; j++, mult *= u) {
          if(j > degree - k) {
            break;
          } else if(j+k < i) {
          } else {
            base += pow(oneOverA, i) * Curve<T, DIM>::comb(degree - i, j + k - i)
                  * pow(-oneOverA, j + k - i) * Curve<T, DIM>::perm(j + k, k)
                  * Curve<T, DIM>::comb(degree, i) * mult;
          }
        }
        basis[i] = base;
      }
      return basis;
    }

    /*
     * Get coefficient matrix of basis functions [B_{0, degree}(u), ..., B_{1, degree}(u)]
     * where the first row of coefficients are for B_{0, degree}(u) and so on.
     *
     * a0 + a1 u + a2 u^2 + ... ad u^d
    */
    Matrix getBasisCoefficientMatrix() const {
      unsigned int degree = m_controlPoints.size() - 1;
      Matrix B(degree+1, degree+1);
      T oneOverA = 1 / m_a;
      for(unsigned int i = 0; i <= degree; i++) {
        // B_{i, degree}(u)
        for(unsigned int j = 0; j <= degree; j++) {
          if(j < i) {
            B(i, j) = 0.0;
          } else {
            B(i, j) = pow(oneOverA, i) * Curve<T, DIM>::comb(degree - i, j - i) *
                      pow(-oneOverA, j - i) * Curve<T, DIM>::comb(degree, i);
          }
        }
      }
      return B;
    }

    inline T parameterSpan() const override {
      return m_a;
    }


    /*
    * Returns the maximum magnitude of kth derivative of the curve
    * with evaluations at every step
    */
    T maxDerivativeMagnitude(unsigned int k, T step) const override {
      T max_mag = std::numeric_limits<T>::lowest();
      for(T t = 0; t < m_a; t += step) {
        max_mag = max(max_mag, eval(t, k).norm());
      }
      return max_mag;
    }
};

}

#endif
