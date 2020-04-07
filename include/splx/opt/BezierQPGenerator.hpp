#ifndef SPLX_BEZIERQPGENERATOR_HPP
#define SPLX_BEZIERQPGENERATOR_HPP
#include <splx/opt/QPGenerator.hpp>
#include <splx/internal/bezier.hpp>
#include <splx/types.hpp>


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

    using Row = typename splx::Row<T>;

    BezierQPGenerator(std::size_t ncpts, T a): Base(ncpts, a) {

    }

    void addIntegratedSquaredDerivativeCost(unsigned int k, T lambda) override {
        if(k >= this->numControlPoints()) return;

        Matrix bern = splx::internal::bezier::bernsteinCoefficientMatrix(this->degree(), this->maxParameter(), k);

        Matrix SQI(this->numControlPoints(), this->numControlPoints());
        SQI.setZero();

        for(Index i = 0; i < this->numControlPoints(); i++) {
            for(Index j = 0; j < this->numControlPoints(); j++) {
                SQI(i, j) = splx::internal::pow(this->maxParameter(), i + j + 1) / (i+j+1);
            }
        }
        
        Matrix cost = 2 * lambda * bern * SQI * bern.transpose();

        Matrix Qcost(this->numControlPoints() * DIM, this->numControlPoints() * DIM);
        Qcost.setZero();

        for(unsigned int i = 0; i < DIM; i++) {
            Qcost.block(i*this->numControlPoints(), 
                        i*this->numControlPoints(), 
                        this->numControlPoints(), 
                        this->numControlPoints()) = cost;
        }

        Base::m_problem.add_Q(Qcost);
    }

    void addEvalCost(T u, unsigned int k, const VectorDIM& target, T lambda) override {
        Row basis = splx::internal::bezier::getBasisRow(this->degree(), this->maxParameter(), u, k);
        Matrix Qext = 2 * lambda * basis.transpose() * basis;
        Matrix Qbig(this->numControlPoints() * DIM, this->numControlPoints() * DIM);
        Qbig.setZero();
        Vector cext = -2 * lambda * basis.transpose();
        Vector cbig(this->numControlPoints() * DIM);
        cbig.setZero();
        for(unsigned int i = 0; i < DIM; i++) {
            Qbig.block(i*this->numControlPoints(), 
                       i*this->numControlPoints(), 
                       this->numControlPoints(), 
                       this->numControlPoints()) = Qext;

            cbig.block(i*this->numControlPoints(), 0, this->numControlPoints(), 1) = cext * target(i);
        }

        Base::m_problem.add_Q(Qbig);
        Base::m_problem.add_c(cbig);
    }

    void addEvalConstraint(T u, unsigned int k, const VectorDIM& target) override {

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