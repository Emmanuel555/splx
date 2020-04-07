#ifndef SPLX_BEZIERQPGENERATOR_HPP
#define SPLX_BEZIERQPGENERATOR_HPP
#include <splx/opt/QPGenerator.hpp>

namespace splx {

template<typename T, unsigned int DIM>
class BezierQPGenerator : public QPGenerator<T, DIM> {
public:
    using Base = QPGenerator<T, DIM>;
    using Index = typename Base::Index;
    using _ParametricCurve = typename Base::_ParametricCurve;
    using Vector = typename Base::Vector;
    using VectorDIM = typename Base::VectorDIM;
    using Hyperplane = typename Base::Hyperplane;
    using Matrix = typename Base::Matrix;

    BezierQPGenerator(std::size_t ncpts, T a): Base(ncpts, a) {

    }

    void addIntegratedSquaredDerivativeCost(unsigned int k, T lambda) override {
        if(k >= this->numControlPoints()) return;
    }

    void addEvalCost(T u, unsigned int k, const VectorDIM& target, T lambda) override {

    }
    void addEvalConstraint(T u, unsigned int k, const VectorDIM& target, T lambda) override {

    }
    void addHyperplaneConstraintAll(const Hyperplane& hp) override {

    }
    void addHyperplaneConstraintAt(T u, const Hyperplane& hp) override {

    }
    void addControlPointLimits(std::size_t i, const VectorDIM& lb, const VectorDIM& ub) override {

    }

public:
    Matrix bernsteinCoefficientMatrixTest(unsigned int k) const { // for test
        return bernsteinCoefficientMatrix(k);
    }

private:

    /*
        Coefficient matrix for the k^th derivative bernstein base functions where each row r
        contains coefficients where k^th derivative of i^th bernstein polynomial of degree d
        can be expressed as r(0) + r(1)u + r(2)u^2 + ... + r(d)u^d
    */
    Matrix bernsteinCoefficientMatrix(unsigned int k) const {
        Matrix bernsteinMtr(this->numControlPoints(), this->numControlPoints());
        bernsteinMtr.setZero();

        if(this->maxParameter() == 0) {
            if(k == 0 && this->numControlPoints() > 0) {
                bernsteinMtr(0, 0) = 1.0;
            }
            return bernsteinMtr;
        }

        unsigned int dcombi = 1;
        for(Index i = 0; 

            i < this->numControlPoints(); 

            dcombi *= (this->degree()-i), 
            dcombi /= (i+1), 
            i++) {

            unsigned int dminicombjmini = 1;
            T min1 = 1;
            T oneOverAPowj = _ParametricCurve::pow(1/this->maxParameter(), i);

            for(Index j = i; 

                j < this->numControlPoints(); 

                dminicombjmini *= (this->degree() - j), 
                dminicombjmini/= (j+1-i), 
                j++, 
                min1 *= -1, 
                oneOverAPowj *= (1/this->maxParameter())) {

                bernsteinMtr(i, j) = dcombi * dminicombjmini * min1 * oneOverAPowj;

            }
        }

        Matrix derivative(this->numControlPoints(), this->numControlPoints());
        derivative.setZero();

        unsigned int jpermk = _ParametricCurve::fac(k);
        for(unsigned int j = k;
            
            j < this->numControlPoints();
            
            jpermk *= (j+1),
            jpermk /= (j+1-k),
            j++) {
            derivative(j, j-k) = jpermk;
        }

        return bernsteinMtr * derivative;
    }
};

} // end namespace splx;

#endif