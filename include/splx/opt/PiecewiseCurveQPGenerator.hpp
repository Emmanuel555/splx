#ifndef SPLX_PIECEWISECURVE_QP_GENERATOR_HPP
#define SPLX_PIECEWISECURVE_QP_GENERATOR_HPP

#include <qp_wrappers/problem.hpp>
#include <splx/opt/QPGenerator.hpp>
#include <splx/opt/BezierQPGenerator.hpp>
#include <splx/types.hpp>
#include <splx/curve/PiecewiseCurve.hpp>
#include <splx/opt/BezierQPGenerator.hpp>
#include <absl/strings/str_cat.h>

namespace splx {

template<typename T, unsigned int DIM>
class PiecewiseCurveQPGenerator {
public:
    using _QPGenerator = QPGenerator<T, DIM>;
    using _BezierQPGenerator = BezierQPGenerator<T, DIM>;
    using _Problem = QPWrappers::Problem<T>;
    using Row = splx::Row<T>;
    using VectorDIM = splx::VectorDIM<T, DIM>;
    using Vector = splx::Vector<T>;
    using Hyperplane = splx::Hyperplane<T, DIM>;
    using AlignedBox = splx::AlignedBox<T, DIM>;
    using Index = splx::Index;
    using Matrix = splx::Matrix<T>;
    using _QPOperations = QPOperations<T, DIM>;
    using _PiecewiseCurve = PiecewiseCurve<T, DIM>;

    PiecewiseCurveQPGenerator() : m_problem(0) {}

    void addPiece(std::shared_ptr<_QPOperations> opt_ptr) {
        m_operations.push_back(opt_ptr);
        
        if(m_cumulativeMaxParameters.empty()) {
            m_cumulativeMaxParameters.push_back(opt_ptr->maxParameter());
        } else {
            m_cumulativeMaxParameters.push_back(
                m_cumulativeMaxParameters.back() + opt_ptr->maxParameter()
            );
        }

        if(m_cumulativeDecisionVars.empty()) {
            m_cumulativeDecisionVars.push_back(opt_ptr->numDecisionVariables());
        } else {
            m_cumulativeDecisionVars.push_back(
                m_cumulativeDecisionVars.back() 
                + opt_ptr->numDecisionVariables()
            );
        }

        m_problem = QPWrappers::Problem<T>(this->numDecisionVariables());
    }

    void setPiece(std::size_t idx, std::shared_ptr<_QPOperations> opt_ptr) {
        m_operations[idx] = opt_ptr;
        this->fixCumulativeStructures(idx);
        m_problem = QPWrappers::Problem<T>(this->numDecisionVariables());
    }

    void removePiece(std::size_t idx) {
        m_operations.erase(m_operations.begin() + idx);
        m_cumulativeDecisionVars.erase(
                m_cumulativeDecisionVars.begin() + idx
        );
        m_cumulativeMaxParameters.erase(
                m_cumulativeMaxParameters.begin() + idx
        );
        this->fixCumulativeStructures(idx);
        m_problem = QPWrappers::Problem<T>(this->numDecisionVariables());
    }

    void removeAllPieces() {
        m_operations.clear();
        m_cumulativeMaxParameters.clear();
        m_cumulativeDecisionVars.clear();
        m_problem = QPWrappers::Problem<T>(0);
    }

    void addBezier(Index ncpts, T a) {
        auto bezptr = std::make_shared<BezierQPGenerator>(ncpts, a);
        auto optptr = std::static_pointer_cast<QPGenerator>(bezptr);
        this->addPiece(optptr);
    }

    void setBezier(std::size_t idx, Index ncpts, T a) {
        auto bezptr = std::make_shared<BezierQPGenerator>(ncpts, a);
        auto optptr = std::static_pointer_cast<QPGenerator>(bezptr);
        this->setPiece(idx, optptr);
    }

    std::size_t numPieces() const {
        return m_cumulativeDecisionVars.size();
    }

    Index numDecisionVariables() const {
        if(m_cumulativeDecisionVars.empty()) return 0;

        return m_cumulativeDecisionVars.back();
    }
    
    void addIntegratedSquaredDerivativeCost(unsigned int k, T lambda) {
        Matrix Q(this->numDecisionVariables(), this->numDecisionVariables());
        Matrix c(this->numDecisionVariables);
        Q.setZero();
        c.setZero();

        for(std::size_t i = 0; i < this->numPieces(); i++) {
            Index dvar_start_idx = (i == 0 ? 0 : m_cumulativeDecisionVars[i-1]);
            Index dvar_count = m_operations[i]->numDecisionVariables();

            auto [Qs, cs] = m_operations[i]->interatedSquaredDerivative(k, lambda);
            
            Q.block(dvar_start_idx, dvar_start_idx, dvar_count, dvar_count) = Qs;
            c.block(dvar_start_idx, 0, dvar_count, 1) = cs;
        }

        m_problem.add_Q(Q);
        m_problem.add_c(c);
    }

    void addEvalCost(T u, unsigned int k, const VectorDIM& target, T lambda) {
        Matrix Q(this->numDecisionVariables(), this->numDecisionVariables());
        Matrix c(this->numDecisionVariables());
        Q.setZero();
        c.setZero();

        auto [idx, param, first_dvar_index, dvar_count] = pieceInfo(u);
        auto [Qs, cs] = m_operations[idx]->evalCost(param, k, target, lambda);
        Q.block(first_dvar_index, first_dvar_index, 
                dvar_count, dvar_count) = Qs;
        c.block(first_dvar_index, 0, dvar_count, 1) = cs;

        m_problem.add_Q(Q);
        m_problem.add_c(c);
    }

    void addEvalConstraint(T u, unsigned int k, const VectorDIM& target) {
        auto [idx, param, first_dvar_index, dvar_count] = pieceInfo(u);
        auto constraints = m_operations[idx]->evalConstraint(param, k, target);

        this->addConstraints(constraints, first_dvar_index, dvar_count);
    }

    void addHyperplaneConstraintForPiece(std::size_t idx, const Hyperplane& hp) {
        Index first_dvar_index = 
                (idx == 0 ? 0 : m_cumulativeDecisionVars[idx-1]);
        Index dvar_count = m_operations[idx]->numDecisionVariables();

        auto constraints = m_operations[idx]->hyperplaneConstraintAll(hp);

        this->addConstraints(constraints, first_dvar_index, dvar_count);
    }

    void addHyperplaneConstraintAll(const Hyperplane& hp) {
        for(std::size_t i = 0; i < this->numPieces(); i++) {
            this->addHyperplaneConstraintForPiece(i, hp);
        }
    }

    void addHyperplaneConstraintAt(T u, const Hyperplane& hp) {
        auto [idx, param, first_dvar_index, dvar_count] = pieceInfo(u);
        auto constraints = m_operations[idx]->hyperplaneConstraintAt(param, hp);
        this->addConstraints(constraints, first_dvar_index, dvar_count);
    }

    void addBoundingBoxConstraint(const AlignedBox& bbox) {
        for(std::size_t i = 0; i < this->numPieces(); i++) {
            auto [lbx, ubx] = m_operations[i]->boundingBoxConstraint(bbox);
            Index first_dvar_index = 
                    (i == 0 ? 0 : m_cumulativeDecisionVars[i-1]);
            Index dvar_count = m_operations[i]->numDecisionVariables();

            assert(lbx.rows() == dvar_count);

            for(Index j = 0; j < dvar_count; j++) {
                m_problem.set_var_limits(first_dvar_index + j, lbx(j), ubx(j));
            }
        }
    }

    // adds continueity constraint between piece idx and piece idx + 1
    // in their kth derivatives
    void addContinuityConstraint(std::size_t idx, unsigned int k) {
        Index first_piece_numdvars = m_operations[idx]->numDecisionVariables();
        Index second_piece_numdvars = 
                m_operations[idx+1]->numDecisionVariables();
        Index first_piece_dvars_start 
            = (idx == 0 ? 0 : m_cumulativeDecisionVars[idx-1]);
        Index second_piece_dvars_start = m_cumulativeDecisionVars[idx];

        for(unsigned int d = 0; d < DIM; d++) {
            Row coeff1 = m_operations[idx]->evalBasisRow(
                            d, m_operations[idx]->maxParameter(), k
            );
            Row coeff2 = m_operations[idx+1]->evalBasisRow(
                            d, m_operations[idx+1]->maxParameter(), k
            );

            Row coeff(this->numDecisionVariables());
            coeff.setZero();
            coeff.block(0, first_piece_dvars_start, 1, first_piece_numdvars) 
                                = coeff1;
            coeff.block(0, second_piece_dvars_start, 1, second_piece_numdvars)
                                = -coeff2;

            m_problem.add_constraint(coeff, 0, 0);
        }
    }

    _PiecewiseCurve extractCurve(const Vector& soln) {
        if(soln.rows() != this->numDecisionVariables()) {
            throw std::domain_error(
                absl::StrCat(
                    "number of decision variables does not match. given: ",
                    soln.rows(),
                    ", required: ",
                    this->numDecisionVariables()
                )
            );
        }

        _PiecewiseCurve piecewise;
        for(std::size_t i = 0; i < this->numPieces(); i++) {
            Index piece_numdvars = m_operations[i]->numDecisionVariables();
            Index piece_dvars_start = 
                    (i==0 ? 0 : m_cumulativeDecisionVars[i-1]);
            
            Vector dvars = soln.block(piece_dvars_start, 0, piece_numdvars, 1);
            piecewise.addPiece(m_operations[i]->extractCurve(dvars));
        }

        return piecewise;
    }

private:
    std::vector<std::shared_ptr<_QPOperations>> m_operations;
    std::vector<T> m_cumulativeMaxParameters;
    std::vector<Index> m_cumulativeDecisionVars;
    QPWrappers::Problem<T> m_problem;

    /*
    * fix cumulative structures starting from the given index
    */
    void fixCumulativeStructures(std::size_t idx) {
        for(std::size_t i = idx; i < this->numPieces(); i++) {
            if(i == 0) {
                m_cumulativeMaxParameters[i] = m_operations[i]->maxParameter();
                m_cumulativeDecisionVars[i] = 
                    m_operations[i]->numDecisionVariables();
            } else {
                m_cumulativeMaxParameters[i] = 
                      m_operations[i]->maxParameter() 
                    + m_cumulativeMaxParameters[i-1];
                m_cumulativeDecisionVars[i] = 
                      m_operations[i]->numDecisionVariables()
                    + m_cumulativeDecisionVars[i-1];
            }
        }
    }

    /*
    * Given a parameter, return the piece index, piece parameter corresponding
    * to the given parameter, first decision variable index and decision
    * variable count for the piece
    */
    std::tuple<std::size_t, T, Index, Index> pieceInfo(T param) const {
        std::size_t idx = std::lower_bound(
                    m_cumulativeMaxParameters.begin(), 
                    m_cumulativeMaxParameters.end(), 
                    param) - m_cumulativeMaxParameters.begin();

        Index first_dvar_index = 0;
        Index dvar_count = m_operations[idx]->numDecisionVariables();
        if(idx != 0) {
            param -= m_cumulativeMaxParameters[idx-1];
            first_dvar_index = m_cumulativeDecisionVars[idx-1];
        }

        return {idx, param, first_dvar_index, dvar_count};
    }

    void addConstraints(const std::vector<std::tuple<Row, T, T>>& cons, 
                        Index first_dvar_index, Index dvar_count) {

        Row coeff(this->numDecisionVariables());

        for(const auto& constraint: cons) {
            coeff.setZero();
            coeff.block(0, first_dvar_index, 1, dvar_count) 
                        = std::get<0>(constraint);
            m_problem.add_constraint(coeff, 
                                     std::get<1>(constraint), 
                                     std::get<2>(constraint)
            );   
        }
    }


};

}

#endif