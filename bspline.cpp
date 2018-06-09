#include "bspline.h"
#include <cassert>
#include <vector>
#include <iostream>
#include <eigen3/Eigen/Dense>

splx::BSpline::BSpline(unsigned int deg, unsigned int dim, double A, double B)
                      : m_degree(deg), m_dimension(dim), m_a(A), m_b(B) {
  assert(m_a < m_b);
}

splx::BSpline::BSpline(unsigned int deg, unsigned int dim, double A, double B,
                 const std::vector<Vec>& cpts)
                : BSpline(deg, dim, A, B) {
  for(int i = 0; i < cpts.size(); i++) {
    assert(cpts[i].rows() == m_dimension);
  }
  m_controlPoints = cpts;
  generateUniformKnotVector();
}

void splx::BSpline::generateUniformKnotVector() {
  assert(m_controlPoints.size() >= m_degree + 1);

  m_knotVector.clear();
  m_knotVector.insert(m_knotVector.begin(), m_degree+1, m_a);
  double insert_count = m_controlPoints.size() - m_degree - 1;
  double step = (m_b - m_a)/(insert_count+1);
  for(int i = 0; i < insert_count; i++)
    m_knotVector.push_back(m_a + (i+1)*step);
  m_knotVector.insert(m_knotVector.end(), m_degree+1, m_b);
}

void splx::BSpline::printKnotVector() const {
  for(int i = 0; i < m_knotVector.size()-1; i++) {
    std::cout << m_knotVector[i] << " ";
  }
  std::cout << m_knotVector[m_knotVector.size() - 1] << std::endl;
}

void splx::BSpline::printControlPoints() const {
  for(int i = 0; i < m_controlPoints.size()-1; i++) {
    std::cout << m_controlPoints[i] << " ";
  }
  std::cout << m_controlPoints[m_controlPoints.size() - 1] << std::endl;
}

void splx::BSpline::printKnotVectorNumbered() const {
  for(int i = 0; i < m_knotVector.size(); i++) {
    std::cout << i << " " << m_knotVector[i] << std::endl;
  }
}

unsigned int splx::BSpline::findSpan(double u) const {
  assert(u >= m_a && u <= m_b);

  if(u == m_b) {
    return m_controlPoints.size() - 1;
  }

  int lo = 0;
  int hi = m_knotVector.size() - 1;
  int mid = (lo + hi) / 2;

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

splx::Vec splx::BSpline::eval(double u) {
  assert(u >= m_a && u <= m_b);

  int js = findSpan(u);

  vector<vector<double> > N(m_controlPoints.size());
  for(int i = 0; i < m_controlPoints.size(); i++) {
    N[i].resize(m_degree+1);
    N[i][0] = (u >= m_knotVector[i] && u < m_knotVector[i+1] ? 1 : 0);
  }

  for(int p = 1; p <= m_degree; p++) {
    for(int j = js - m_degree; j <= js; j++) {
      N[j][p] = 
    }
  }

}

splx::Vec splx::BSpline::eval(double u, int n) {
  assert(u >= m_a && u <= m_b);
  if(n == 0) {
    return eval(u);
  }
}
