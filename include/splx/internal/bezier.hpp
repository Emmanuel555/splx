#ifndef SPLX_INTERNAL_BEZIER_H
#define SPLX_INTERNAL_BEZIER_H
#include <Eigen/Dense>
#include <splx/types.hpp>
#include <splx/internal/combinatorics.hpp>

namespace splx {
namespace internal {
namespace bezier {

/*
* given a bezier curve degree and max parameter, return row r so that
* multiplying the row with control points gives the kth derivative
* of the curve at u
*
* f^k(u) = \sum_{i=0}^degree r(i) * p(i)
* return r
*/
template<typename T>
Row<T> getBasisRow(unsigned int degree, T maxParameter, T u, unsigned int k) {
    if(u < 0 || u > maxParameter) {
    throw std::domain_error(
        std::string("u is outside of the range [0, ")
        + std::to_string(maxParameter)
        + std::string("]")
    );
    }

    if(maxParameter == 0) {
        Row<T> result(degree+1);
        result.setZero();
        if(k == 0 && degree >= 0) {
            result(0) = 1.0;
        }
        return result;
    }

    Row<T> result(degree + 1);

    T oneOverA = 1/maxParameter;
    for(unsigned int i = 0; i <= degree; i++) {
    T base = 0.0;
    T mult = 1.0;
    for(unsigned int j = 0; j+k <= degree; j++, mult *= u) {
        if(j+k >= i) {
        // base += pow(oneOverA, i) * this->comb(degree-i, j+k-i)
        //       * pow(-oneOverA, j+k-i) * this->perm(j+k, k)
        //       * mult;

        base += splx::internal::comb(degree-i, j+k-i)
            * splx::internal::pow(oneOverA, j+k) * splx::internal::perm(j+k, k)
            * mult * ((j+k-i)%2 == 0 ? 1 : -1);
        }
    }
    base *= splx::internal::comb(degree, i);
    result(i) = base;
    }
    return result;
}

/*
    Coefficient matrix for the k^th derivative bernstein base functions where each row r
    contains coefficients where k^th derivative of i^th bernstein polynomial of degree d
    can be expressed as r(0) + r(1)u + r(2)u^2 + ... + r(d)u^d
*/
template<typename T>
Matrix<T> bernsteinCoefficientMatrix(unsigned int degree, T maxParameter, unsigned int k) {
    Matrix<T> bernsteinMtr(degree+1, degree+1);
    bernsteinMtr.setZero();

    if(maxParameter == 0) {
        if(k == 0 && degree >= 0) {
            bernsteinMtr(0, 0) = 1.0;
        }
        return bernsteinMtr;
    }

    unsigned int dcombi = 1;
    for(Index i = 0; 

        i < degree+1; 

        dcombi *= (degree-i), 
        dcombi /= (i+1), 
        i++) {

        unsigned int dminicombjmini = 1;
        T min1 = 1;
        T oneOverAPowj = splx::internal::pow(1/maxParameter, i);

        for(Index j = i; 

            j < degree+1; 

            dminicombjmini *= (degree - j), 
            dminicombjmini/= (j+1-i), 
            j++, 
            min1 *= -1, 
            oneOverAPowj *= (1/maxParameter)) {

            bernsteinMtr(i, j) = dcombi * dminicombjmini * min1 * oneOverAPowj;

        }
    }

    Matrix<T> derivative(degree+1, degree+1);
    derivative.setZero();

    unsigned int jpermk = splx::internal::fac(k);
    for(unsigned int j = k;
        
        j < degree+1;
        
        jpermk *= (j+1),
        jpermk /= (j+1-k),
        j++) {
        derivative(j, j-k) = jpermk;
    }

    return bernsteinMtr * derivative;
}

} // namespace bezier
} // namespace internal
} // namespace splx

#endif