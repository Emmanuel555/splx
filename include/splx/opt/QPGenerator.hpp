#ifndef SPLX_QPGENERATOR_HPP
#define SPLX_QPGENERATOR_HPP
#include <qp_wrappers/problem.hpp>
#include <splx/curve/ParametricCurve.hpp>

namespace splx {
    
template<typename T, unsigned int DIM>
class QPGenerator {
public:
    using _ParametricCurve = ParametricCurve<T, DIM>;
    using VectorDIM = typename _ParametricCurve::VectorDIM;
    using Hyperplane = typename _ParametricCurve::Hyperplane;

    using Problem = QPWrappers::Problem<T>;
    using Index = typename Problem::Index;

    using Matrix = typename _ParametricCurve::Matrix;
    using Vector = typename _ParametricCurve::Vector;

    QPGenerator(std::size_t n, T a): m_ncpts(n), m_a(a), m_problem(n * DIM) {

    }

    /**
    *   adds the cost lambda * \int_{0}^{m_a} ||df^k(u)/du^k||_2^2 du
    */
    virtual void addIntegratedSquaredDerivativeCost(unsigned int k, T lambda) = 0;

    /**
    *   adds the cost lambda * ||f^k(u)-target||^2
    */
    virtual void addEvalCost(T u, unsigned int k, const VectorDIM& target, T lambda) = 0;

    /**
    *   adds constraint that f^k(u) = target
    */
    virtual void addEvalConstraint(T u, unsigned int k, const VectorDIM& target) = 0;

    /*
    * given hyperplane hp(x) = a_1x_1+a_2x_2+...+a_nx_n + d = 0 
    * such that n=(a_1, ..., a_n) is the normal,
    * add a constraint so that hp(f(u)) <= 0 for all u 
    */
    virtual void addHyperplaneConstraintAll(const Hyperplane& hp) = 0;

    /*
    * given hyperplane hp(x) = a_1x_1+a_2x_2+...+a_nx_n + d = 0 
    * such that n=(a_1, ..., a_n) is the normal,
    * add a constraint so that hp(f(u)) <= 0 at u
    */
    virtual void addHyperplaneConstraintAt(T u, const Hyperplane& hp) = 0;
    virtual void addControlPointLimits(std::size_t i, const VectorDIM& lb, const VectorDIM& ub) = 0;

    std::size_t numControlPoints() const {
        return m_ncpts;
    }
    
    unsigned int degree() const {
        return m_ncpts - 1;
    }

    T maxParameter() const {
        return m_a;
    }

    const Problem& getProblem() const {
        return m_problem;
    }

protected:
    std::size_t m_ncpts;
    T m_a;
    QPWrappers::Problem<T> m_problem;
}; // end QPGenerator


} // end namespace splx


#endif