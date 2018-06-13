#include "bspline.h"
#include <cassert>
#include <vector>
#include <iostream>
#include <eigen3/Eigen/Dense>

using std::cout;
using std::endl;

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
