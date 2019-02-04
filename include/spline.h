#ifndef SPLX_SPLINE_H
#define SPLX_SPLINE_H


#include <vector>
#include <memory>
#include "curve.h"
#include <numeric>
#include "bezier.h"
#include "bspline.h"
#include <iostream>

using std::cout;
using std::endl;
using std::vector;

namespace splx {

template<typename T, unsigned int DIM>
class Spline {
  public:
    using VectorDIM = typename Eigen::Matrix<T, DIM, 1>;
    using Vector = typename Eigen::Matrix<T, Eigen::Dynamic, 1>;
    using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using Hyperplane = Eigen::Hyperplane<T, DIM>;
    using MatrixDIM = Eigen::Matrix<T, DIM, DIM, Eigen::RowMajor>;
    using AlignedBox = Eigen::AlignedBox<T, DIM>;

    using QPMatrices = typename Curve<T, DIM>::QPMatrices;

    class OutOfBoundsException : public std::runtime_error {
      public:
        OutOfBoundsException(const char *msg): std::runtime_error(msg) {

        }
    };

    /*
     *
    */
    std::vector<QPMatrices> getQPMatrices() const {
      std::vector<QPMatrices> result;
      for(unsigned int i = 0; i < m_pieces.size(); i++) {
        result.push_back(m_pieces[i]->getQPMatrices());
      }
      return result;
    }


    /*
     * Evaluate k^th derivative of spline at parameter u
    */
    VectorDIM eval(T u, unsigned int k) const {
      assert(u <= totalSpan() && u >= 0);
      std::pair<unsigned int, T> info = curveInfo(u);
      return m_pieces[info.first]->eval(info.second, k);
    }

    /*
     * Returns which bezier given u corresponds to with the parameter
     * on the bezier.
    */
    std::pair<unsigned int, T> curveInfo(T u) const {
      unsigned int i = 0;
      while(i < m_pieces.size() && u > m_pieces[i]->parameterSpan()) {
        u -= m_pieces[i]->parameterSpan();
        i++;
      }

      if(i == m_pieces.size()) {
         i = m_pieces.size() - 1;
         u = m_pieces[i]->parameterSpan();
         // throw OutOfBoundsException("given parameter is out of bounds");
      }

      return std::make_pair(i, u);
    }

    T totalSpan() const {
      return std::accumulate(m_pieces.begin(), m_pieces.end(), T(),
                            [] (T acc, const std::shared_ptr<Curve<T, DIM> >& crv) -> T {
                              return acc + crv->parameterSpan();
                            }
             );
    }

    void addPiece(const std::shared_ptr<Curve<T, DIM> >& piece) {
      m_pieces.push_back(piece);
    }

    void addPiece(const BSpline<T, DIM>& piece) {
      std::shared_ptr<Curve<T, DIM> > ptr = std::make_shared<BSpline<T, DIM> >(piece);
      m_pieces.push_back(ptr);
    }

    void addPiece(const Bezier<T, DIM>& piece) {
      std::shared_ptr<Curve<T, DIM> > ptr = std::make_shared<Bezier<T, DIM> >(piece);
      m_pieces.push_back(ptr);
    }

    void addPiece(const Curve<T, DIM>& piece) {
      std::shared_ptr<Curve<T, DIM> > ptr = std::make_shared<Curve<T, DIM> >(piece);
      m_pieces.push_back(ptr);
    }

    std::shared_ptr<Curve<T, DIM> > getPiece(unsigned int i) const {
      return m_pieces[i];
    }

    unsigned int numPieces() const {
      return m_pieces.size();
    }

    void extendQPIntegratedSquaredDerivative(std::vector<QPMatrices>& QPs, unsigned int k, T lambda) const {
      assert(QPs.size() == m_pieces.size());

      for(unsigned int i = 0; i < m_pieces.size(); i++) {
        m_pieces[i]->extendQPIntegratedSquaredDerivative(QPs[i], k, lambda);
      }
    }

    void extendQPPositionAt(std::vector<QPMatrices>& QPs, T u, const VectorDIM& pos, T theta) const {
      assert(QPs.size() == m_pieces.size());

      std::pair<unsigned int, T> info = curveInfo(u);
      m_pieces[info.first]->extendQPPositionAt(QPs[info.first], info.second, pos, theta);
    }

    void extendQPBeginningConstraint(std::vector<QPMatrices>& QPs, unsigned int k, const VectorDIM& target) const {
      assert(QPs.size() == m_pieces.size());

      m_pieces[0]->extendQPBeginningConstraint(QPs[0], k, target);
    }

    void extendQPDecisionConstraint(std::vector<QPMatrices>& QPs, T lb, T ub) const {
      assert(QPs.size() == m_pieces.size());

      for(unsigned int i = 0; i < QPs.size(); i++) {
        m_pieces[i]->extendQPDecisionConstraint(QPs[i], lb, ub);
      }
    }

    /*
     * require i^th curve to be in the negative side of hyperplane hp
    */
    void extendQPHyperplaneConstraint(std::vector<QPMatrices>& QPs, unsigned int i, const Hyperplane& hp) const {
      assert(QPs.size() == m_pieces.size());
      assert(i < QPs.size());

      m_pieces[i]->extendQPHyperplaneConstraint(QPs[i], hp);
    }

    QPMatrices combineQPMatrices(const std::vector<QPMatrices>& QPs) {
      QPMatrices QP;
      unsigned int varCount = 0;
      for(unsigned int i = 0; i < QPs.size(); i++) {
        varCount += QPs[i].x.rows();
      }

      QP.H.resize(varCount, varCount);
      QP.g.resize(varCount);
      QP.lbX.resize(varCount);
      QP.ubX.resize(varCount);
      QP.x.resize(varCount);



      for(unsigned int i = 0; i < varCount; i++) {
        for(unsigned int j = 0; j < varCount; j++) {
          QP.H(i, j) = 0.0;
        }
        QP.g(i) = 0.0;
        QP.lbX(i) = 0.0;
        QP.ubX(i) = 0.0;
        QP.x(i) = 0.0;
      }

      QP.A.resize(0, varCount);
      QP.lbA.resize(0);
      QP.ubA.resize(0);

      unsigned int lastIdx = 0;
      for(unsigned int i = 0; i < QPs.size(); i++) {
        unsigned int varC = QPs[i].x.rows();
        QP.H.block(lastIdx, lastIdx, varC, varC) = QPs[i].H;
        QP.g.block(lastIdx, 0, varC, 1) = QPs[i].g;
        QP.lbX.block(lastIdx, 0, varC, 1) = QPs[i].lbX;
        QP.ubX.block(lastIdx, 0, varC, 1) = QPs[i].ubX;
        QP.x.block(lastIdx, 0, varC, 1) = QPs[i].x;

        unsigned int ridx = QP.A.rows();

        QP.A.conservativeResize(QP.A.rows() + QPs[i].A.rows(), varCount);
        QP.lbA.conservativeResize(QP.lbA.rows() + QPs[i].lbA.rows());
        QP.ubA.conservativeResize(QP.ubA.rows() + QPs[i].ubA.rows());

        for(unsigned int j = 0; j < QPs[i].A.rows(); j++) {
          QP.lbA(ridx+j) = QPs[i].lbA(j);
          QP.ubA(ridx+j) = QPs[i].ubA(j);
          for(unsigned int k = 0; k < varCount; k++) {
            QP.A(ridx + j, k) = 0.0;
          }
          QP.A.block(ridx+j, lastIdx, 1, QPs[i].x.rows()) = QPs[i].A.block(j, 0, 1, QPs[i].x.rows());
        }
        lastIdx += varC;
      }

      return QP;
    }

    /*
     * This must be called after combine operation
     * QP is the combined matrix, QPs are individual matrices returned from getQPMatrices
     *
     * Imposes k^th degree continuity constraint between piece i and i+1
    */
    void extendQPContinuityConstraint(QPMatrices& QP, const std::vector<QPMatrices>& QPs, unsigned int i, unsigned int k) {
      assert(i < QPs.size() - 1);
      const unsigned int varCount = QP.x.rows();
      unsigned int varCountBeforeI = 0;
      for(unsigned int j = 0; j < i; j++) {
        varCountBeforeI += QPs[j].x.rows();
      }

      typename Curve<T, DIM>::Row iRow = m_pieces[i]->getQPBasisRow(m_pieces[i]->parameterSpan(), k);
      typename Curve<T, DIM>::Row ipRow = m_pieces[i+1]->getQPBasisRow(0, k);
      ipRow = -1 * ipRow;

      unsigned int ridx = QP.A.rows();

      QP.A.conservativeResize(QP.A.rows() + DIM, QP.A.cols());
      QP.lbA.conservativeResize(QP.lbA.rows() + DIM);
      QP.ubA.conservativeResize(QP.ubA.rows() + DIM);

      for(unsigned int d = 0; d < DIM; d++) {
        QP.lbA(ridx + d) = 0.0;
        QP.ubA(ridx + d) = 0.0;
        for(unsigned int j = 0; j < varCount; j++) {
          QP.A(ridx+d, j) = 0.0;
        }
        QP.A.block(ridx+d, varCountBeforeI + d * iRow.cols(), 1, iRow.cols()) = iRow;
        QP.A.block(ridx+d, varCountBeforeI + QPs[i].x.rows() + d * ipRow.cols(), 1, ipRow.cols()) = ipRow;
      }

    }

    /*
     * QP is the qp matrices of combined optimization
     * Gets the real data from QP matrix
     * Fixes the QPs matrices in the way as well
     * which can be nice, dunno
    */
    void loadControlPoints(const QPMatrices& QP, std::vector<QPMatrices>& QPs) {
      assert(QPs.size() == m_pieces.size());
      unsigned int start = 0;
      for(unsigned int i = 0; i < QPs.size(); i++) {
        unsigned int varcount = QPs[i].x.rows();
        QPs[i].x = QP.x.block(start, 0, varcount, 1);
        m_pieces[i]->loadControlPoints(QPs[i]);
        start += varcount;
      }

    }

    /*
    * Returns true if the spline is inside of the aligned box at any point from
    * parameter from to parameter to
    */
    bool intersects(const AlignedBox& box, T from, T to, T radius) const {
      const T step = 0.01;
      VectorDIM rad(radius, radius, radius);
      for(T t = from; t <= to; t += step) {
        auto info = curveInfo(t);
        unsigned int i = info.first;
        T u = info.second;
        VectorDIM pos = m_pieces[i]->eval(u, 0);
        AlignedBox robot_box(pos - rad, pos + rad);
        if(box.intersects(robot_box)) {
          return true;
        }
      }
      return false;
    }

    /*
    * Returns true if the spline is inside of any of the aligned boxes at any point
    * from parameter from to parameter to
    */
    bool intersects(const vector<AlignedBox>& boxes, T from, T to, T radius) const {
      for(const auto& box: boxes) {
        if(intersects(box, from, to, radius)) {
          return true;
        }
      }
      return false;
    }

    /*
    * Returns true if the ith piece is in the negative side of hp
    */
    bool onNegativeSide(const Hyperplane& hp, unsigned int i) const {
      return m_pieces[i]->onNegativeSide(hp);
    }

    /*
    * Returns true if the spline is the negative side of the hyperplane
    */
    bool onNegativeSide(const Hyperplane& hp) const {
      for(int i = 0; i < m_pieces.size(); i++) {
        if(!onNegativeSide(hp, i))
          return false;
      }
      return true;
    }

    /*
    * Returns the max magnitude of k^th derivative along the curve
    * Using evaluations at every step
    */
    T maxDerivativeMagnitude(unsigned int k, T step) const {
      T max_mag = std::numeric_limits<T>::lowest();
      for(const auto& piece: m_pieces) {
        max_mag = max(max_mag, piece->maxDerivativeMagnitude(k, step));
      }

      return max_mag;
    }


  private:
    std::vector<std::shared_ptr<Curve<T, DIM> > > m_pieces;


};

}


#endif
