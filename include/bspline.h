#ifndef SPLX_BSPLINE_H
#define SPLX_BSPLINE_H
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#include <vector>
#include "spline.h"
#include <utility>

namespace splx {

template<typename T, unsigned int DIM>
class BSpline : public Spline<T, DIM> {
public:

  using Vector = typename Spline<T, DIM>::Vector;
  using Hyperplane = typename Spline<T, DIM>::Hyperplane;
  using Matrix = typename Spline<T, DIM>::Matrix;
  using VectorDIM = typename Spline<T, DIM>::VectorDIM;
  using QPMatrices = typename Spline<T, DIM>::QPMatrices;

  /**
   * Construct a b-spline where basis functions are of
   * degree deg and curve has parameter u in [A, B]
   *
   * @fails if A > B
  */
  BSpline(unsigned int deg, T A, T B): m_degree(deg), m_a(A), m_b(B) {
    assert(m_a <= m_b);
  }

  /**
   * Construct a b-spline where basis functions are of
   * degree deg and curve has parameter u \in [A, B] with given initial
   * control points
   *
   * @fails if A >= B
   * @fails if the size of m_controlPoints is less than m_degree+1
  */
  BSpline(unsigned int deg, T A,
          T B, const std::vector<VectorDIM>& cpts) : m_controlPoints(cpts),
          m_degree(deg), m_a(A), m_b(B) {
    assert(m_a <= m_b);
    generateClampedUniformKnotVector();
  }


  /**
    * Evaluates the k^{th} derivative of the spline at u \in [0,1]
    *
    * @fails if u is not in [m_a, m_b]
  */
  VectorDIM eval(T u, unsigned int k) const override {
    assert(u >= m_a && u <= m_b);
    VectorDIM result;
    for(unsigned int i = 0; i < DIM; i++) {
      result(i) = 0.0;
    }

    unsigned int je = findSpan(u);
    unsigned int js = je < m_degree ? 0 : je - m_degree;

    std::vector<T> N = evalBasisFuncs(u, m_degree, k, js, je);

    for(unsigned int j = js; j <= je; j++) {
      result += m_controlPoints[j] * N[j - js];
    }

    return result;
  }


  /**
   * In all functions that extends these matrices, order of variables is assumed to be
   * p0x, p1x, p2x, ..., pnx, p0y, p1y, ..., pny, ...
   *
  */
  QPMatrices getQPMatrices() const override {
    QPMatrices QP;
    unsigned int S = m_controlPoints.size() * DIM;
    QP.H.resize(S, S);
    QP.g.resize(S);
    for(unsigned int i = 0; i < S; i++) {
      QP.g(i) = 0.0;
      for(unsigned int j = 0; j < S; j++)
        QP.H(i, j) = 0.0;
    }
    QP.A.resize(0, S);
    QP.lbA.resize(0);
    QP.ubA.resize(0);

    QP.lbX.resize(S);
    QP.ubX.resize(S);
    for(unsigned int i = 0; i < S; i++) {
      QP.lbX(i) = std::numeric_limits<T>::lowest();
      QP.ubX(i) = std::numeric_limits<T>::max();
    }

    QP.x.resize(S);
    for(unsigned int i = 0; i < m_controlPoints.size(); i++) {
      for(unsigned int d = 0; d < DIM; d++) {
        QP.x(d*m_controlPoints.size() + i) = m_controlPoints[i](d);
      }
    }
    return QP;
  }


  /**
   * Add integral from m_a to m_b of square of norm of k^th derivative of the spline
   * to the hessian matrix H with scalar lambda
   *
   * order of variables is assumed to be
   * p0x, p1x, p2x, ..., pnx, p0y, p1y, ..., pny, ...
  */
  void extendQPIntegratedSquaredDerivative(QPMatrices& QP, unsigned int k, T lambda) const {
    if(k > m_degree)
      return;

    Matrix D(m_degree+1, m_degree+1); // get k^th derivative coefficients from a0 + a1u + ... a(m_degree)u^m_degree
    for(unsigned int m = 0; m < m_degree+1; m++) {
      for(unsigned int n = 0; n < m_degree+1; n++) {
        D(m, n) = 0.0;
      }
      if(m <= m_degree-k) {
        D(m, m+k) = Spline<T, DIM>::perm(m+k, k);
      }
    }


    for(unsigned int j = 0; j < m_knotVector.size() - 1; j++) {
      // integrate from m_knotVector[j] to m_knotVector[j+1]
      if(m_knotVector[j] == m_knotVector[j+1])
        continue;
      unsigned int js = j < m_degree ? 0U : j - m_degree;
      Matrix M = getBasisCoefficientMatrix(js, j, m_degree, j).transpose();
      Matrix Mext(m_degree+1, m_controlPoints.size());
      for(unsigned int m = 0; m <= m_degree; m++) {
        for(unsigned int n = 0; n < m_controlPoints.size(); n++) {
          if(n >= js && n <= j) {
            Mext(m, n) = M(m, n-js);
          } else {
            Mext(m, n) = 0.0;
          }
        }
      }
      Matrix SQI(m_degree+1, m_degree+1); // get the integral of the square of the polynomial.
      for(unsigned m = 0; m <= m_degree; m++) {
        for(unsigned int n = 0; n <= m_degree; n++) {
          SQI(m, n) = 2.0 * (std::pow(m_knotVector[j+1], m+n+1) - std::pow(m_knotVector[j], m+n+1)) / (m+n+1);
        }
      }

      Matrix Hext = lambda * Mext.transpose() * D.transpose() * SQI * D * Mext;
      for(unsigned int d = 0; d < DIM; d++) {
        QP.H.block(d*m_controlPoints.size(), d*m_controlPoints.size(), m_controlPoints.size(), m_controlPoints.size()) += Hext;
      }
    }
  }

  /**
   * Add the cost theta * ||f(u) - pos||^2 to H and g.
   *
   * order of variables is assumed to be
   * p0x, p1x, p2x, ..., pnx, p0y, p1y, ..., pny, ...
  */

  void extendQPPositionAt(QPMatrices& QP, T u, const VectorDIM& pos, T theta) const {
    assert(u >= m_a && u <= m_b);
    unsigned int je = findSpan(u);
    unsigned int js = je < m_degree ? 0 : je - m_degree;

    std::vector<T> basis = evalBasisFuncs(u, m_degree, 0, js, je);
    Vector Mext(m_controlPoints.size());
    for(unsigned int i = 0; i<m_controlPoints.size(); i++) {
      if(i >= js && i<=je) {
        Mext(i) = basis[i-js];
      } else {
        Mext(i) = 0.0;
      }
    }

    Matrix Hext = 2 * theta * Mext * Mext.transpose();

    for(unsigned int d = 0; d < DIM; d++) {
      Vector Gext = -2 * theta * pos[d] * Mext;
      QP.g.block(d*m_controlPoints.size(), 0, m_controlPoints.size(), 1) += Gext;
      QP.H.block(d*m_controlPoints.size(), d*m_controlPoints.size(), m_controlPoints.size(), m_controlPoints.size()) += Hext;
    }
  }

  /**
   * Add constraint that requires the k^th derivative of spline at u=0 to be target.
  */
  void extendQPBeginningConstraint(QPMatrices& QP, unsigned int k, const VectorDIM& target) const {
    assert(k <= m_degree);
    unsigned int ridx = QP.A.rows();
    QP.A.conservativeResize(QP.A.rows() + DIM, QP.A.cols());
    QP.lbA.conservativeResize(QP.lbA.rows() + DIM);
    QP.ubA.conservativeResize(QP.ubA.rows() + DIM);

    unsigned int je = findSpan(m_a);
    unsigned int js = je < m_degree ? 0 : je - m_degree;

    std::vector<T> basis = evalBasisFuncs(m_a, m_degree, k, js, je);
  /*
    for(int i = 0; i < basis.size(); i++) {
      cout << basis[i] << " " << endl;
    }
    int a; cin >> a;
  */
    unsigned int S = m_controlPoints.size() * DIM;


    for(unsigned int d = 0; d < DIM; d++) {
      for(unsigned int m = 0; m < S; m++) {
        if(m >= d*m_controlPoints.size()+js && m <= d*m_controlPoints.size()+je) {
          QP.A(ridx+d, m) = basis[m - d*m_controlPoints.size() - js];
        } else {
          QP.A(ridx+d, m) = 0.0;
        }
      }
      QP.lbA(ridx+d) = target(d);
      QP.ubA(ridx+d) = target(d);
    }
  }

  /**
   * Add constraints that requires control points from index from to index to to be on the
   * negative side of hp.
  */
  void extendQPHyperplaneConstraint(QPMatrices& QP, unsigned int from, unsigned int to, const Hyperplane& hp) const {
    assert(to < m_controlPoints.size());

    unsigned int ridx = QP.A.rows();
    QP.A.conservativeResize(QP.A.rows() + to-from+1, QP.A.cols());
    QP.lbA.conservativeResize(QP.lbA.rows() + to-from+1);
    QP.ubA.conservativeResize(QP.ubA.rows() + to-from+1);


    unsigned int S = DIM * m_controlPoints.size();

    for(unsigned int i = from; i<=to; i++) {
      for(unsigned int s = 0; s < S; s++) {
        QP.A(ridx + i - from, s) = 0.0;
      }
      for(unsigned int d = 0; d < DIM; d++) {
        QP.A(ridx + i - from, d*m_controlPoints.size() + i) = hp.normal()(d);
      }
      QP.lbA(ridx + i - from) = std::numeric_limits<T>::lowest();
      QP.ubA(ridx + i - from) = -hp.offset();
    }
  }

  /**
   * Add constraint that requires curve to be on the negative side of hp
   * when from <= u <= to
   *
   * Effectively, enforces all points that effects the curve in [from, to] to be
   * in the negative side of the hp.
  */
  void extendQPHyperplaneConstraint(QPMatrices& QP, T from, T to, const Hyperplane& hp) const {
    assert(from >= m_a && from <= m_b);
    assert(to >= m_a && to <= m_b);

    auto aff = affectingPoints(from, to);

    extendQPHyperplaneConstraint(QP, aff.first, aff.second, hp);
  }

  /**
   * Add constraints that requires all decision variables to be less than ub and more than lb.
   *
   * Effectively, it constraints the curve to be inside a box.
  */
  void extendQPDecisionConstraint(QPMatrices&QP, T lb, T ub) const {
    for(unsigned int i=0; i < QP.lbX.rows(); i++) {
      QP.lbX(i) = lb;
      QP.ubX(i) = ub;
    }
  }
  /*
   * returns the inclusive index range of points that effects the curve
   * from u = from to u = to.
  */
  std::pair<unsigned int, unsigned int> affectingPoints(T from, T to) const {
    assert(to >= from && to <= m_b && from >= m_a);
    unsigned int js = findSpan(from);
    unsigned int je = findSpan(to);

    js = js < m_degree ? 0 : js - m_degree;

    return std::make_pair(js, je);
  }


  /*
   * Getter for k^th control point
  */
  const VectorDIM& getCP(unsigned int k) const {
    return m_controlPoints[k];
  }

  /*
   * Interpolate from point 'from' to point 'to' with n points.
   * 'from' is the 1st point, 'to' is nth point.
   * Repeat 'to' m_degree+1 times to make the curve end at 'to'.
   * In the end there are n + m_degree points added.
   *
   * returns inclusive index range of added points in m_controlPoints array.
  */
  std::pair<unsigned int, unsigned int> interpolateEndAtTo(const VectorDIM& from, const VectorDIM& to, unsigned int n) {
    assert(n >= 2);
    std::pair<unsigned int, unsigned int> result;
    result.first = m_controlPoints.size();
    VectorDIM step = (to - from) / (n-1);
    for(unsigned int i = 0; i < n; i++) {
      m_controlPoints.push_back(from + step * i);
    }

    for(unsigned int i = 0; i < m_degree; i++) {
      m_controlPoints.push_back(to);
    }
    result.second = m_controlPoints.size() - 1;
    return result;
  }


  /**
    * Generates uniform knot vector from scratch that is clamped
    *
    * @assumes control points are already set
  */
  void generateClampedUniformKnotVector() {
    assert(m_controlPoints.size() >= m_degree + 1);
    m_knotVector.clear();
    m_knotVector.insert(m_knotVector.begin(), m_degree+1, m_a);
    unsigned int insert_count = m_controlPoints.size() - m_degree - 1;
    T step = (m_b - m_a)/(insert_count+1);
    for(int i = 0; i < insert_count; i++)
      m_knotVector.push_back(m_a + (i+1)*step);
    m_knotVector.insert(m_knotVector.end(), m_degree+1, m_b);
  }

  /**
    * Generates uniform knot vector from scratch that is not clamped
    *
    * @assumes control points are already set
    *
  */
  void generateNonclampedUniformKnotVector() {
    m_knotVector.clear();
    unsigned int insert_count = m_controlPoints.size() + m_degree + 1;
    T step = (m_b - m_a) / (insert_count - 1);
    for(unsigned int i = 0; i < insert_count - 1; i++) {
      m_knotVector.push_back(m_a + i * step);
    }
    m_knotVector.push_back(m_b); // seperated for numerical reasons
  }


  /**
   * Generates non uniform non clamped knot vector from scratch
   * according to given weights w.
   *
   * divites total span to |w| number of intervals
   * where the lengths of each interval i is proportional to w[i]
   * and each interval has the same number of knots
  */
  void generateNonclampedNonuniformKnotVector(const std::vector<T>& w) {
    m_knotVector.clear();
    unsigned int insert_count = m_controlPoints.size() + m_degree + 1;
    T span = m_b - m_a;
    T step = span / (insert_count - 1);
    m_knotVector.push_back(m_a); // seperated for numerical issues
    for(unsigned int i = 1; i < insert_count - 1; i++) {
      m_knotVector.push_back(m_a + step * i);
    }
    m_knotVector.push_back(m_b); // seperated for numerical issues
  }

  /**
    * Generates non uniform knot vector from scratch
    * according to given weights w.
    *
    * divides total span to |w| number of intervals
    * where the length of each interval i is proportional to w[i]
    * and each interval has the same number of knots.
    */
  void generateClampedNonuniformKnotVector(const std::vector<T>& w) {
    assert(m_controlPoints.size() >= m_degree + 1);
    m_knotVector.clear();
    m_knotVector.insert(m_knotVector.begin(), m_degree + 1, m_a);
    unsigned int insert_count = m_controlPoints.size() - m_degree;
    unsigned int insert_per_step = insert_count / w.size();
    T totalw = std::accumulate(w.begin(), w.end(), 0.0);
    T last_end = m_a;
    for(size_t i = 0; i < w.size() - 1; i++) {
      T ratio = w[i] / totalw;
      T cur_end = last_end + (m_b - m_a) * ratio;
      T step = (cur_end - last_end) / insert_per_step;
      for(unsigned int j = 0; j < insert_per_step; j++) {
        m_knotVector.push_back(last_end + step * (j+1));
      }
      last_end = cur_end;
    }

    /* last element is handled seperately to deal with under/overflows */
    insert_per_step = insert_count - insert_per_step * (w.size() - 1);
    T step = (m_b - last_end) / insert_per_step;
    for(unsigned int j = 0; j < insert_per_step; j++) {
      m_knotVector.push_back(last_end + step * (j+1));
    }

    m_knotVector.insert(m_knotVector.end(), m_degree, m_b);
  }

  /**
   * Clear control points array.
  */
  void clearControlPoints() {
    m_controlPoints.clear();
  }

  /**
   * Load control points from QP.x
  */
  void loadControlPoints(const QPMatrices& QP) {
    for(unsigned int i = 0; i < m_controlPoints.size(); i++) {
      for(unsigned int d = 0; d < DIM; d++) {
        m_controlPoints[i](d) = QP.x(d * m_controlPoints.size() + i);
      }
    }
  }

  /**
    DBG FUNCTIONS
  */
  VectorDIM eval_dbg(T u) const {
    unsigned int je = findSpan(u);
    unsigned int js = je < m_degree ? 0 : je - m_degree;
    Matrix mtr = getBasisCoefficientMatrix(js, je, m_degree, je);

    Vector uvec(m_degree+1);
    T uy = 1;
    for(unsigned int i=0; i<m_degree+1; i++) {
      uvec(i) = uy;
      uy *= u;
    }

    Vector basis = mtr * uvec;

    VectorDIM res;
    for(unsigned int i = 0;i < DIM;i++)
      res(i) = 0.0;

    for(unsigned int i = 0; i < m_degree+1; i++) {
      res += basis(i) * m_controlPoints[i+js];
    }
    return res;
  }
  void printKnotVector() const {
    for(size_t i = 0; i < m_knotVector.size()-1; i++) {
      std::cout << m_knotVector[i] << " ";
    }
    std::cout << m_knotVector[m_knotVector.size() - 1] << std::endl;
  }

  void printControlPoints() const {
    for(size_t i = 0; i < m_controlPoints.size()-1; i++) {
      std::cout << m_controlPoints[i] << " ";
    }
    std::cout << m_controlPoints[m_controlPoints.size() - 1] << std::endl;
  }
  void printKnotVectorNumbered() const {
    for(size_t i = 0; i < m_knotVector.size(); i++) {
      std::cout << i << " " << m_knotVector[i] << std::endl;
    }
  }


  unsigned int m_degree; // degree of basis functions
  T m_a; // first p+1 knot values
  T m_b; // last p+1 knot values
  std::vector<VectorDIM> m_controlPoints; // control points

private:
  /**
    dimension that spline is defined in,
    defines the dimension of control points as well.
  */
  std::vector<T> m_knotVector; // knot vector


  /**
    * Finds the index i of m_knotVector where u falls into [u_i, u_{i+1})
    *
    * @fails if u not in [0,1]
  */
  unsigned int findSpan(T u) const {
    assert(u >= m_a && u <= m_b);

    if(u == m_b) {
      unsigned int idx = m_knotVector.size() - 1;
      while(m_knotVector[idx] == u) {
        idx--;
      }
      return idx;
    }
    auto it = std::upper_bound(m_knotVector.begin(), m_knotVector.end(), u);
    return it - m_knotVector.begin() - 1;
  }

  /**
   * Evaluate k^th derivative of basis functions [N_{from,deg}, ..., N_{to,deg}] at u.
  */
  std::vector<T> evalBasisFuncs(T u, unsigned int deg, unsigned int k, unsigned int from, unsigned int to) const {
    //cout << "evaluating for u = " << u << endl;
    if(from > to) {
      return std::vector<T>();
    }

    if(k > deg) {
      return std::vector<T>(to - from + 1, 0.0);
    }

    std::vector<std::vector<T> > N(2);
    N[0].resize(to + deg - from + 1, 0.0);
    N[1].resize(to + deg - from + 1, 0.0);

    unsigned int end_span = findSpan(m_b);

    for(unsigned int j = from; j <= to + deg; j++) {
      N[0][j - from] = (u >= m_knotVector[j] && u < m_knotVector[j+1] ? 1.0 : 0.0);
      if(u == m_b && j == end_span) {
        N[0][j - from] = 1.0;
      }
      //cout << N[0][j - from] << endl;
    }
    //int a; cin >> a;

    for(unsigned int p = 1; p <= deg - k; p++) {
      int i = p & 0x1;
      int pi = (p - 1) & 0x1;
      for(unsigned int j = from; j <= to + deg - p; j++) {
        N[i][j-from] =
          (N[pi][j-from] == 0.0 ? 0.0 : N[pi][j-from] * (u - m_knotVector[j]) / (m_knotVector[j+p] - m_knotVector[j]))
        + (N[pi][j-from+1] == 0.0 ? 0.0 : N[pi][j-from+1] * (m_knotVector[j+p+1] - u) / (m_knotVector[j+p+1] - m_knotVector[j+1]));
        //cout << N[i][j-from] << endl;
      }
      //cin >> a;
    }

    for(unsigned int  p = deg - k + 1; p <= deg; p++) {
      int i = p & 0x1;
      int pi = (p - 1) & 0x1;
      for(unsigned int j = from; j <= to  + deg - p; j++) {
        N[i][j-from] = p * ((N[pi][j-from] == 0.0 ? 0.0 : N[pi][j-from] / (m_knotVector[j+p] - m_knotVector[j]))
                     - (N[pi][j-from+1] == 0.0 ? 0.0 : N[pi][j-from+1] / (m_knotVector[j+p+1] - m_knotVector[j+1])));
      }
    }

    //N[deg & 0x1].resize(to - from + 1);
    return N[deg & 0x1];
  }

  /**
   * Get coefficient matrix of basis functions [N_{from,p}(u) ... N_{to,p}(u)] in interval
   * [m_knotVector[i], m_knotVector[i+1]) where first row is the coefficients of N_{from, p}
   * and last row is the coefficients of N_{to, p}
   *
   * a0 + a1u + a2u^2 + ... + apu^p
  */
  Matrix getBasisCoefficientMatrix(unsigned int from, unsigned int to, unsigned int p, unsigned int i) const {
    Matrix result(to-from+1, p+1);
    for(unsigned int ix = 0; ix < to-from+1; ix++) {
      for(unsigned int iy = 0; iy < p+1; iy++) {
        result(ix, iy) = 0.0;
      }
    }

    if(p == 0) {
      for(unsigned int j = from; j <= to; j++) {
        result(j-from, 0) = (j == i ? 1.0 : 0.0);
      }
    } else {
      Matrix prevpower = getBasisCoefficientMatrix(from ,to+1, p-1, i);

      for(unsigned int j = from; j<=to; j++) {
        result(j-from, 0) += (m_knotVector[j+p] == m_knotVector[j] ? 0.0 :
            (-m_knotVector[j] * prevpower(j-from, 0) / (m_knotVector[j+p] - m_knotVector[j])));
        result(j-from, 0) += (m_knotVector[j+p+1] == m_knotVector[j+1] ? 0.0 :
            (m_knotVector[j+p+1] * prevpower(j-from+1, 0) / (m_knotVector[j+p+1] - m_knotVector[j+1])));

        result(j-from, p) += (m_knotVector[j+p] == m_knotVector[j] ? 0.0 :
            (prevpower(j-from, p-1) / (m_knotVector[j+p] - m_knotVector[j])));
        result(j-from, p) += (m_knotVector[j+p+1] == m_knotVector[j+1] ? 0.0 :
            (-prevpower(j-from+1, p-1) / (m_knotVector[j+p+1] - m_knotVector[j+1])));
        for(unsigned int k = 1; k < p; k++) {
          result(j-from, k) += (m_knotVector[j+p] == m_knotVector[j] ? 0.0:
                          (prevpower(j-from, k-1) / (m_knotVector[j+p] - m_knotVector[j])));
          result(j-from, k) += (m_knotVector[j+p] == m_knotVector[j] ? 0.0:
                          (-m_knotVector[j] * prevpower(j-from, k) / (m_knotVector[j+p] - m_knotVector[j])));
          result(j-from, k) += (m_knotVector[j+p+1] == m_knotVector[j+1] ? 0.0:
                          (-prevpower(j+1-from, k-1) / (m_knotVector[j+p+1] - m_knotVector[j+1])));
          result(j-from, k) += (m_knotVector[j+p+1] == m_knotVector[j+1] ? 0.0:
                          (m_knotVector[j+p+1] * prevpower(j+1-from, k) / (m_knotVector[j+p+1] - m_knotVector[j+1])));
        }
      }
    }

    return result;
  }

};

}

#endif
