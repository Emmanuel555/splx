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
    using ControlPoints = typename Base::ControlPoints;

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
        Row basis = splx::internal::bezier::getBasisRow(this->degree(), this->maxParameter(), u, k);
        for(unsigned int i = 0; i < DIM; i++) {
            Row coeff(this->numControlPoints() * DIM);
            coeff.setZero();
            coeff.block(0, i*this->numControlPoints(), 1, this->numControlPoints()) = basis;
            Base::m_problem.add_constraint(coeff, target(i), target(i));
        }
    }

    void addHyperplaneConstraintAll(const Hyperplane& hp) override {
        for(Index i = 0; i < Base::m_ncpts; i++) {
            Row coeff(this->numControlPoints() * DIM);
            coeff.setZero();
            for(Index j = 0; j < DIM; j++) {
                coeff(j * this->numControlPoints() + i) = hp.normal()(j);
            }
            Base::m_problem.add_constraint(coeff, std::numeric_limits<T>::lowest(), -hp.offset());
        }
    }

    void addHyperplaneConstraintAt(T u, const Hyperplane& hp) override {
        Row basis = splx::internal::bezier::getBasisRow(this->degree(), this->maxParameter(), u, 0);

        Row coeff(this->numControlPoints() * DIM);
        for(Index i = 0; i < DIM; i++) {
            coeff.block(0, i*this->numControlPoints(), 1, this->numControlPoints()) = basis * hp.normal()(i);
        }
        Base::m_problem.add_constraint(coeff, std::numeric_limits<T>::lowest(), -hp.offset());
    }

    void addControlPointLimits(const VectorDIM& lb, const VectorDIM& ub) override {
        for(Index j = 0; j < DIM; j++) {
            for(Index i = 0; i < this->numControlPoints(); i++) {
                Base::m_problem.set_var_limits(j * this->numControlPoints() + i, lb(j), ub(j));
            }
        }
    }

    ControlPoints controlPoints(const Vector& solution) const override {
        ControlPoints cpts;
        for(unsigned int i = 0; i < Base::m_ncpts; i++) {
            Vector cpt(DIM);
            for(unsigned int j = 0; j < DIM; j++) {
                cpt(j) = solution(j * this->numControlPoints() + i);
            }
        }
    }

private:

    unsigned int degree() const {
        return this->numControlPoints() - 1;
    }

};

} // end namespace splx;

#endif