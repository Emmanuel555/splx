#include "spline.h"
#include <algorithm>

unsigned int splx::Spline::fac(unsigned int n) const {
  unsigned int res = 1;
  for(unsigned int i = 2; i<=n; i++)
    res *= i;
  return res;
}


unsigned int splx::Spline::comb(unsigned int n, unsigned int k) const {
  k = std::min(k, n-k);
  unsigned int top = 1;
  unsigned int bottom = 1;
  for(unsigned int i=0; i<k; i++) {
    bottom *= (i+1);
    top *= (n-i);
  }
  return top / bottom;
}


unsigned int splx::Spline::perm(unsigned int n, unsigned int k) const {
  return comb(n, k) * fac(k);
}


splx::Matrix splx::Spline::convertHessianToUpperTriangular(const splx::Matrix& H) const {
  unsigned int S = H.rows();
  splx::Matrix U(H);
  for(unsigned int i = 0; i < S; i++) {
    U(i, i) /= 2.0;
    for(unsigned int j = 0; j < i; j++) {
      U(i, j) = 0.0;
    }
  }
  return U;
}


void splx::Spline::convertHessianToUpperTriangular(splx::Matrix& H) const {
  unsigned int S = H.rows();
  for(unsigned int i = 0; i < S; i++) {
    H(i, i) /= 2.0;
    for(unsigned int j = 0; j < i; j++) {
      H(i, j) = 0.0;
    }
  }
}
