#ifndef SPLX_QPGENERATOR_HPP
#define SPLX_QPGENERATOR_HPP
#include <qp_wrappers/problem.hpp>
#include <splx/curve/ParametricCurve.hpp>

namespace splx {
    
template<typename T, unsigned int DIM>
class QPGenerator {
public:
    using ParametricCurve = ParametricCurve<T, DIM>;
    using VectorDIM = ParametricCurve::VectorDIM;
    using Hyperplane = ParametricCurve::Hyperplane;

    using Problem = QPWrappers::Problem<T>;
    using Index = Problem::Index;

    QPGenerator(std::size_t n): ncpts(n), problem(n * DIM) {

    }

    virtual void addIntegratedSquaredDerivativeCost(unsigned int k, T lambda) = 0;
    virtual void addEvalCost(T u, unsigned int k, const VectorDIM& target, T lambda) = 0;
    virtual void addEvalConstraint(T u, unsigned int k, const VectorDIM& target, T lambda) = 0;
    virtual void addHyperplaneConstraintAll(const Hyperplane& hp) = 0;
    virtual void addHyperplaneConstraintAt(T u, const Hyperplane& hp) = 0;
    virtual void addControlPointLimits(std::size_t i, const VectorDIM& lb, const VectorDIM& ub) = 0;

    std::size_t numControlPoints() const {
        return ncpts;
    }
    
protected:
    std::size_t ncpts;
    QPWrappers::Problem<T> problem;
}; // end QPGenerator


} // end namespace splx


#endif