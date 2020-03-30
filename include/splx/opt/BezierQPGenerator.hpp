#ifndef SPLX_BEZIERQPGENERATOR_HPP
#define SPLX_BEZIERQPGENERATOR_HPP
#include <splx/opt/QPGenerator.hpp>

namespace splx {

template<typename T, unsigned int DIM>
class BezierQPGenerator : public QPGenerator<T, DIM> {
public:
    using Base = QPGenerator<T, DIM>;
    BezierQPGenerator(std::size_t ncpts): Base(ncpts) {

    }

    void addIntegratedSquaredDerivativeCost(unsigned int k, T lambda) override {

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


private:

};

} // end namespace splx;

#endif