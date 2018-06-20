#include "bspline.h"
#include <cassert>
#include <vector>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <cmath>

using std::cout;
using std::endl;
using std::cin;

splx::BSpline::BSpline(unsigned int deg, unsigned int dim, double A, double B)
                      : m_degree(deg), m_dimension(dim), m_a(A), m_b(B) {
  assert(m_a <= m_b);
}

splx::BSpline::BSpline(unsigned int deg, unsigned int dim, double A, double B,
                 const std::vector<Vec>& cpts)
                : BSpline(deg, dim, A, B) {
  for(size_t i = 0; i < cpts.size(); i++) {
    assert(cpts[i].rows() == m_dimension);
  }
  m_controlPoints = cpts;
  assert(m_controlPoints.size() >= m_degree + 1);
  generateUniformKnotVector();
}

void splx::BSpline::generateUniformKnotVector() {
  m_knotVector.clear();
  m_knotVector.insert(m_knotVector.begin(), m_degree+1, m_a);
  double insert_count = m_controlPoints.size() - m_degree - 1;
  double step = (m_b - m_a)/(insert_count+1);
  for(int i = 0; i < insert_count; i++)
    m_knotVector.push_back(m_a + (i+1)*step);
  m_knotVector.insert(m_knotVector.end(), m_degree+1, m_b);
}

void splx::BSpline::printKnotVector() const {
  for(size_t i = 0; i < m_knotVector.size()-1; i++) {
    std::cout << m_knotVector[i] << " ";
  }
  std::cout << m_knotVector[m_knotVector.size() - 1] << std::endl;
}

void splx::BSpline::printControlPoints() const {
  for(size_t i = 0; i < m_controlPoints.size()-1; i++) {
    std::cout << m_controlPoints[i] << " ";
  }
  std::cout << m_controlPoints[m_controlPoints.size() - 1] << std::endl;
}

void splx::BSpline::printKnotVectorNumbered() const {
  for(size_t i = 0; i < m_knotVector.size(); i++) {
    std::cout << i << " " << m_knotVector[i] << std::endl;
  }
}

unsigned int splx::BSpline::findSpan(double u) const {
  assert(u >= m_a && u <= m_b);

  if(u == m_b) { // special cases
    return m_controlPoints.size() - 1;
  }

  unsigned int lo = 0;
  unsigned int hi = m_knotVector.size() - 1;
  unsigned int mid = (lo + hi) / 2;

  while(u < m_knotVector[mid] || u >= m_knotVector[mid+1]) {
    if(u < m_knotVector[mid]) {
      hi = mid;
    } else {
      lo = mid+1;
    }
    mid = (hi + lo) / 2;
  }
  return mid;
}

std::vector<double> splx::BSpline::evalBasisFuncs(double u, unsigned int deg, unsigned int k, unsigned int from, unsigned int to) const {
  if(from > to) {
    return std::vector<double>();
  }

  if(k > deg) {
    return std::vector<double>(to - from + 1, 0.0);
  }

  std::vector<std::vector<double> > N(2);
  N[0].resize(to + deg - from + 1, 0.0);
  N[1].resize(to + deg - from + 1, 0.0);


  for(unsigned int j = from; j <= to + deg; j++) {
    N[0][j - from] = (u >= m_knotVector[j] && u < m_knotVector[j+1] ? 1 : 0);
    if(u == m_b && j == m_controlPoints.size() - 1) {
      N[0][j - from] = 1.0;
    }
  }

  for(unsigned int p = 1; p <= deg - k; p++) {
    int i = p & 0x1;
    int pi = (p - 1) & 0x1;
    for(unsigned int j = from; j <= to + deg - p; j++) {
      N[i][j-from] =
        (N[pi][j-from] == 0.0 ? 0.0 : N[pi][j-from] * (u - m_knotVector[j]) / (m_knotVector[j+p] - m_knotVector[j]))
      + (N[pi][j-from+1] == 0.0 ? 0.0 : N[pi][j-from+1] * (m_knotVector[j+p+1] - u) / (m_knotVector[j+p+1] - m_knotVector[j+1]));
    }
  }

  for(unsigned int  p = deg - k + 1; p <= deg; p++) {
    int i = p & 0x1;
    int pi = (p - 1) & 0x1;
    for(unsigned int j = from; j <= to  + deg - p; j++) {
      N[i][j-from] = p * ((N[pi][j-from] == 0.0 ? 0.0 : N[pi][j-from] / (m_knotVector[j+p] - m_knotVector[j]))
                   - (N[pi][j-from+1] == 0.0 ? 0.0 : N[pi][j-from+1] / (m_knotVector[j+p+1] - m_knotVector[j+1])));
    }
  }

  N[deg & 0x1].resize(to - from + 1);
  return N[deg & 0x1];
}


splx::Vec splx::BSpline::eval(double u, unsigned int k) const {
  assert(u >= m_a && u <= m_b);
  Vec result(m_dimension);
  for(unsigned int i = 0; i < m_dimension; i++) {
    result(i) = 0.0;
  }

  unsigned int je = findSpan(u);

  std::vector<double> N = evalBasisFuncs(u, m_degree, k, je - m_degree, je);

  for(unsigned int j = je - m_degree; j <= je; j++) {
    result += m_controlPoints[j] * N[j - je + m_degree];
  }

  return result;
}

splx::QPMatrices splx::BSpline::getQPMatrices() const {
  QPMatrices QP;
  unsigned int S = m_controlPoints.size() * m_dimension;
  QP.H.resize(S, S);
  QP.g.resize(S);
  for(unsigned int i = 0; i < S; i++) {
    QP.g(i) = 0.0;
    for(unsigned int j = 0; j < S; j++)
      QP.H(i, j) = 0.0;
  }
  QP.A.resize(0, S);
  QP.lb.resize(0);
  QP.ub.resize(0);

  return QP;
}



void splx::BSpline::extendQPIntegratedSquaredDerivative(QPMatrices& QP, unsigned int k, double lambda) const {
  if(k > m_degree)
    return;

  Matrix D(m_degree+1, m_degree+1); // get k^th derivative coefficients from a0 + a1u + ... a(m_degree)u^m_degree
  for(unsigned int m = 0; m < m_degree+1; m++) {
    for(unsigned int n = 0; n < m_degree+1; n++) {
      D(m, n) = 0.0;
    }
    if(m <= m_degree-k) {
      D(m, m+k) = perm(m+k, k);
    }
  }


  for(unsigned int j = m_degree; j < m_controlPoints.size(); j++) {
    // integrate from m_knotVector[j] to m_knotVector[j+1]
    Matrix M = getBasisCoefficientMatrix(j - m_degree, j, m_degree, j).transpose();
    Matrix Mext(m_degree+1, m_controlPoints.size());
    for(unsigned int m = 0; m <= m_degree; m++) {
      for(unsigned int n = 0; n < m_controlPoints.size(); n++) {
        if(n >= j - m_degree && n <= j) {
          Mext(m, n) = M(m, n-j+m_degree);
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
    for(unsigned int d = 0; d < m_dimension; d++) {
      QP.H.block(d*m_controlPoints.size(), d*m_controlPoints.size(), m_controlPoints.size(), m_controlPoints.size()) += Hext;
    }
  }
}


void splx::BSpline::extendQPPositionAt(QPMatrices& QP, double u, const splx::Vec& pos, double theta) const {
  assert(u >= m_a && u <= m_b);
  unsigned int je = findSpan(u);

  std::vector<double> res = evalBasisFuncs(u, m_degree, 0, je-m_degree, je);
  Vec Mext(m_controlPoints.size());
  for(unsigned int i = 0; i<m_controlPoints.size(); i++) {
    if(i >= je-m_degree && i<=je) {
      Mext(i) = res[i-je+m_degree];
    } else {
      Mext(i) = 0.0;
    }
  }

  Matrix Hext = 2 * theta * Mext * Mext.transpose();

  for(unsigned int d = 0; d < m_dimension; d++) {
    Vec Gext = -2 * pos[d] * Mext;
    QP.g.block(d*m_controlPoints.size(), 0, m_controlPoints.size(), 1) += Gext;
    QP.H.block(d*m_controlPoints.size(), d*m_controlPoints.size(), m_controlPoints.size(), m_controlPoints.size()) += Hext;
  }
}


splx::Matrix splx::BSpline::getBasisCoefficientMatrix(unsigned int from, unsigned int to, unsigned int p, unsigned int i) const {
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


splx::Vec splx::BSpline::eval_dbg(double u) const {
  unsigned int je = findSpan(u);
  Matrix mtr = getBasisCoefficientMatrix(je-m_degree, je, m_degree, je);

  Vec uvec(m_degree+1);
  double uy = 1;
  for(unsigned int i=0; i<m_degree+1; i++) {
    uvec(i) = uy;
    uy *= u;
  }

  Vec basis = mtr * uvec;

  Vec res(m_dimension);
  for(unsigned int i=0;i<m_dimension;i++)
    res(i) = 0.0;

  for(unsigned int i = 0; i<m_degree+1; i++) {
    res += basis(i) * m_controlPoints[i+je-m_degree];
  }
  return res;
}
