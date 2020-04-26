#ifndef SPLX_INTERNAL_COMBINATORICS_H
#define SPLX_INTERNAL_COMBINATORICS_H

#include <cmath>
#include <algorithm>

namespace splx {
namespace internal {

unsigned int fac(unsigned int n) {
    unsigned int res = 1;
    for(unsigned int i = 2; i<=n; i++)
        res *= i;
    return res;
}

unsigned int comb(unsigned int n, unsigned int k) {
    k = std::min(k, n-k);
    unsigned int top = 1;
    unsigned int bottom = 1;
    for(unsigned int i=0; i<k; i++) {
        bottom *= (i+1);
        top *= (n-i);
    }
    return top / bottom;
}

unsigned int perm(unsigned int n, unsigned int k) {
    return fac(n) / fac(n-k);
}

template<typename T>
T pow(T base, unsigned int exp) {
    if(base == 0 && exp == 0) {
        return 1;
    }
    return std::pow(base, exp);
}

} // namespace internal
} // namespace splx

#endif