#include "spline.h"
#include <algorithm>

unsigned int splx::Spline::fac(unsigned int n) const {
  unsigned int res = 1;
  for(int i = 2; i<=n; i++)
    res *= i;
  return res;
}


unsigned int splx::Spline::comb(unsigned int n, unsigned int k) const {
  k = std::min(k, n-k);
  int top = 1;
  int bottom = 1;
  for(int i=0; i<k; i++) {
    bottom *= (i+1);
    top *= (n-i);
  }
  return top / bottom;
}


unsigned int splx::Spline::perm(unsigned int n, unsigned int k) const {
  return comb(n, k) * fac(k);
}
