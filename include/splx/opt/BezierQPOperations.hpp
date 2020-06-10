#ifndef SPLX_BEZIERQPGENERATOR_HPP
#define SPLX_BEZIERQPGENERATOR_HPP
#include <splx/opt/QPOperations.hpp>
#include <splx/internal/bezier.hpp>
#include <splx/types.hpp>
#include <absl/strings/str_cat.h>
#include <splx/curve/Bezier.hpp>


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
    using _ParametricCurve = ParametricCurve<T, DIM>;
    using _Bezier = Bezier<T, DIM>;
    using Constraint = splx::Constraint<T>;

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


    std::vector<Constraint> evalConstraint(
                T u, unsigned int k, const VectorDIM& target,
                bool soft_convertible = false,
                T soft_convertible_weight = T(1)) const override {

        _Base::parameterCheck(u);

        Row basis = splx::internal::bezier::getBasisRow(
                                    this->degree(), 
                                    _Base::maxParameter(), 
                                    u, 
                                    k
        );

        std::vector<Constraint> constraints;
        for(unsigned int i = 0; i < DIM; i++) {
            Row coeff(_Base::numDecisionVariables());
            coeff.setZero();
            coeff.block(0, i*this->numControlPoints(), 
                        1, this->numControlPoints()) = basis;
            constraints.emplace_back(coeff, target(i), target(i),
                    soft_convertible, soft_convertible_weight);
        }
        return constraints;
    }

    std::vector<Constraint> hyperplaneConstraintAll(
                    const Hyperplane& hp,
                    bool soft_convertible = false,
                    T soft_convertible_weight = T(1)) const override {

        std::vector<Constraint> constraints;
        for(Index i = 0; i < this->numControlPoints(); i++) {
            Row coeff(_Base::numDecisionVariables());
            coeff.setZero();
            for(Index j = 0; j < DIM; j++) {
                coeff(j * this->numControlPoints() + i) = hp.normal()(j);
            }
            constraints.emplace_back(coeff, 
                                     std::numeric_limits<T>::lowest(), 
                                     -hp.offset(),
                                     soft_convertible, soft_convertible_weight
            );
        }
        return constraints;
    };

    std::vector<Constraint> hyperplaneConstraintAt(
                                    T u, const Hyperplane& hp,
                                    bool soft_convertible = false,
                                    T soft_convertible_weight = T(1))
                                    const override {

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

        return {{coeff, std::numeric_limits<T>::lowest(), -hp.offset(),
                        soft_convertible, soft_convertible_weight}};
    }

    std::pair<Vector, Vector> boundingBoxConstraint(
                        const AlignedBox& bbox) const override {
        Vector lbx(_Base::numDecisionVariables());
        Vector ubx(_Base::numDecisionVariables());
        for(Index j = 0; j < DIM; j++) {
            for(Index i = 0; i < this->numControlPoints(); i++) {
                lbx(j * this->numControlPoints() + i) = bbox.min()(j);
                ubx(j * this->numControlPoints() + i) = bbox.max()(j);
            }
        }
        return std::make_pair(lbx, ubx);
    };

    Vector getDVarsForSegment(
            const VectorDIM& from, const VectorDIM& to) const override {
        VectorDIM step = (to - from) / (this->numControlPoints()-1);
        Vector dvars(this->numDecisionVariables());
        for(Index i = 0; i < this->numControlPoints(); i++) {
            VectorDIM cpt = from + step * i;
            if(i == 0) cpt = from;
            if(i == this->numControlPoints() - 1) cpt = to;

            for(unsigned int j = 0; j < DIM; j++) {
                dvars(j *  this->numControlPoints() + i) = cpt(j);
            }
        }
        return dvars;
    }

    Index numControlPoints() const { 
        return m_ncpts; 
    }
    void numControlPoints(Index ncpts) {
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

    std::shared_ptr<_ParametricCurve> extractCurve(
            const Vector& soln) const override { 

        if(soln.rows() != _Base::numDecisionVariables()) {
            throw std::domain_error(
                absl::StrCat(
                    "number of decision variables does not match. given: ",
                    soln.rows(),
                    ", required: ",
                    _Base::numDecisionVariables()
                )
            );
        }

        typename _Bezier::ControlPoints cpts;
        for(Index i = 0; i < this->numControlPoints(); i++) {
            VectorDIM cpt;
            for(Index j = 0; j < DIM; j++) {
                cpt(j) = soln(j*this->numControlPoints() + i);
            }
            cpts.push_back(cpt);
        }

        auto bezptr = std::make_shared<_Bezier>(_Base::maxParameter(), cpts);

        return std::static_pointer_cast<_ParametricCurve>(bezptr);
    }

private:
    Index m_ncpts;
};


} // end namespace splx;

#endif