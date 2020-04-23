#ifndef SPLX_BEZIERQPGENERATOR_HPP
#define SPLX_BEZIERQPGENERATOR_HPP
#include <splx/opt/QPGenerator.hpp>
#include <splx/internal/bezier.hpp>
#include <splx/types.hpp>


namespace splx {

template<typename T, unsigned int DIM>
class BezierQPOperations : public QPOperations<T, DIM> {
public:
    using _Base = QPOperations<T, DIM>;
    using VectorDIM = splx::VectorDIM<T, DIM>;
    using Row = splx::Row<T>;
    using Vector = splx::Vector<T>;
    using Index = splx::Index;
    using Matrix = splx::Matrix<T>;
    using Hyperplane = splx::Hyperplane<T, DIM>;
    using AlignedBox = splx::AlignedBox<T, DIM>;

    BezierQPOperations(Index ncpts, T a) 
        : _Base(ncpts * DIM, a), m_ncpts(ncpts) {

    }

    ~BezierQPOperations() {

    }

    Row evalBasisRow(unsigned int d, T u, unsigned int k) const override {
        _Base::parameterCheck(u);
        _Base::dimensionIndexCheck(d);

        Row coeff(_Base::numDecisionVariables());
        coeff.setZero();

        Row r = splx::internal::bezier::getBasisRow(
                    this->degree(), 
                    _Base::maxParameter(), 
                    u, 
                    k
        );

        coeff.block(0, d*this->numControlPoints(), 
                    1, this->numControlPoints()) = r;

        return coeff;
    };

    std::pair<Matrix, Vector> 
    integratedSquaredDerivativeCost(unsigned int k, T lambda) const override {

        Matrix Q(_Base::numDecisionVariables(), _Base::numDecisionVariables());
        Q.setZero();
        Vector c(_Base::numDecisionVariables());
        c.setZero();


        if(k <= this->degree()) {
            Matrix bern = splx::internal::bezier::bernsteinCoefficientMatrix(
                                    this->degree(), 
                                    _Base::maxParameter(), 
                                    k
            );

            Matrix SQI(this->numControlPoints(), this->numControlPoints());
            SQI.setZero();

            for(Index i = 0; i < this->numControlPoints(); i++) {
                for(Index j = 0; j < this->numControlPoints(); j++) {
                    SQI(i, j) = splx::internal::pow(
                                    _Base::maxParameter(), 
                                    i + j + 1) / (i+j+1);
                }
            }
            
            Matrix cost = 2 * lambda * bern * SQI * bern.transpose();

            

            for(unsigned int i = 0; i < DIM; i++) {
                Q.block(i*this->numControlPoints(), 
                            i*this->numControlPoints(), 
                            this->numControlPoints(), 
                            this->numControlPoints()) = cost;
            }

        }

        
        return std::make_pair(Q, c);
    };


    std::pair<Matrix, Vector> evalCost(T u, unsigned int k, 
                                       const VectorDIM& target, 
                                       T lambda) const override {

        _Base::parameterCheck(u);

        Row basis = splx::internal::bezier::getBasisRow(
                            this->degree(), 
                            _Base::maxParameter(), 
                            u, 
                            k
        );

        Matrix Qext = 2 * lambda * basis.transpose() * basis;
        Matrix Qbig(_Base::numDecisionVariables(), 
                    _Base::numDecisionVariables());

        Qbig.setZero();
        Vector cext = -2 * lambda * basis.transpose();
        Vector cbig(_Base::numDecisionVariables());
        cbig.setZero();

        for(unsigned int i = 0; i < DIM; i++) {
            Qbig.block(i*this->numControlPoints(), 
                       i*this->numControlPoints(), 
                       this->numControlPoints(), 
                       this->numControlPoints()) = Qext;

            cbig.block(i*this->numControlPoints(), 
                       0, this->numControlPoints(), 1) = cext * target(i);
        }

        return std::make_pair(Qbig, cbig);
    }


    std::vector<std::tuple<Row, T, T>> evalConstraint(
                T u, unsigned int k, const VectorDIM& target) const override {

        _Base::parameterCheck(u);

        Row basis = splx::internal::bezier::getBasisRow(
                                    this->degree(), 
                                    _Base::maxParameter(), 
                                    u, 
                                    k
        );

        std::vector<std::tuple<Row, T, T>> constraints;
        for(unsigned int i = 0; i < DIM; i++) {
            Row coeff(_Base::numDecisionVariables());
            coeff.setZero();
            coeff.block(0, i*this->numControlPoints(), 
                        1, this->numControlPoints()) = basis;
            constraints.emplace_back(coeff, target(i), target(i));
        }
        return constraints;
    }

    std::vector<std::tuple<Row, T, T>> hyperplaneConstraintAll(
                    const Hyperplane& hp) const override {

        std::vector<std::tuple<Row, T, T>> constraints;
        for(Index i = 0; i < this->numControlPoints(); i++) {
            Row coeff(_Base::numDecisionVariables());
            coeff.setZero();
            for(Index j = 0; j < DIM; j++) {
                coeff(j * this->numControlPoints() + i) = hp.normal()(j);
            }
            constraints.emplace_back(coeff, 
                                     std::numeric_limits<T>::lowest(), 
                                     -hp.offset()
            );
        }
        return constraints;
    };

    std::vector<std::tuple<Row, T, T>> hyperplaneConstraintAt(
                                    T u, const Hyperplane& hp) const override {

        _Base::parameterCheck(u);

        Row basis = splx::internal::bezier::getBasisRow(
                            this->degree(), 
                            _Base::maxParameter(), 
                            u, 
                            0
        );

        Row coeff(this->numControlPoints() * DIM);
        for(Index i = 0; i < DIM; i++) {
            coeff.block(0, i*this->numControlPoints(), 
                        1, this->numControlPoints()) = basis * hp.normal()(i);
        }

        return {{coeff, std::numeric_limits<T>::lowest(), -hp.offset()}};
    }

    std::pair<Vector, Vector> boundingBoxConstraint(
                        const AlignedBox& bbox) const override {
        Vector lbx(_Base::numDecisionVariables());
        Vector ubx(_Base::numDecisionVariables());
        for(Index j = 0; j < DIM; j++) {
            for(Index i = 0; i < this->numControlPoints(); i++) {
                lbx(j * this->numControlPoints() + i, bbox.min()(j));
                ubx(j * this->numControlPoints() + i, bbox.max()(j));
            }
        }
        return std::make_pair(lbx, ubx);
    };

    Index numControlPoints() const { 
        return m_ncpts; 
    }
    void numControlPoints(Index ncpts) const { 
        m_ncpts = ncpts; 
        _Base::numDecisionVariables(ncpts * DIM); 
    }
    Index degree() const {
        if(this->numControlPoints() == 0) {
            throw std::domain_error(
                absl::StrCat(
                    "degree not defined for a bezier "
                     "curve with 0 control points"
                )
            );
        }
        
        return this->numControlPoints() - 1;
    }

private:
    Index m_ncpts;
};

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
    using AlignedBox = splx::AlignedBox<T, DIM>;

    using _BezierQPOperations = BezierQPOperations<T, DIM>;

    BezierQPGenerator(Index ncpts, T a): 
                    Base(ncpts * DIM), 
                    m_operations(ncpts, a) {

    }


    void addIntegratedSquaredDerivativeCost(unsigned int k, T lambda) override {
        auto [Q, c] = m_operations.integratedSquaredDerivativeCost(k, lambda);
        Base::m_problem.add_Q(Q);
        Base::m_problem.add_c(c);
    }

    void addEvalCost(T u, unsigned int k, 
                     const VectorDIM& target, T lambda) override {
        auto [Q, c] = m_operations.evalCost(u, k, target, lambda);
        Base::m_problem.add_Q(Q);
        Base::m_problem.add_c(c);
    }

    void addEvalConstraint(T u, unsigned int k, 
                           const VectorDIM& target) override {
        auto constraints = m_operations.evalConstraint(u, k, target);
        Base::addConstraints(constraints);
    }

    void addHyperplaneConstraintAll(const Hyperplane& hp) override {
        auto constraints = m_operations.hyperplaneConstraintAll(hp);
        Base::addConstraints(constraints);
    }

    void addHyperplaneConstraintAt(T u, const Hyperplane& hp) override {
        auto constraints = m_operations.hyperplaneConstraintAt(u, hp);
        Base::addConstraints(constraints);
    }

    void addBoundingBoxConstraint(const AlignedBox& bbox) override {
        auto [lbx, ubx] = m_operations.boundingBoxConstraint(bbox);
        for(Index i = 0; i < m_operations.numDecisionVariables(); i++) {
            Base::m_problem.set_var_limits(i, lbx(i), ubx(i));
        }
    }

private:
    _BezierQPOperations m_operations;

};

} // end namespace splx;

#endif