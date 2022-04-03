#ifndef SPLX_PIECEWISECURVE_QP_GENERATOR_HPP
#define SPLX_PIECEWISECURVE_QP_GENERATOR_HPP

#include <qp_wrappers/problem.hpp>
#include <splx/opt/QPOperations.hpp>
#include <splx/opt/BezierQPOperations.hpp>
#include <splx/types.hpp>
#include <splx/curve/PiecewiseCurve.hpp>
#include <splx/opt/BezierQPOperations.hpp>
#include <absl/strings/str_cat.h>
#include <Eigen/StdVector>

namespace splx {

template<typename T, unsigned int DIM>
class PiecewiseCurveQPGenerator {
public:
    using _BezierQPOperations = BezierQPOperations<T, DIM>;
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
    using StdVectorVectorDIM
        = std::vector<VectorDIM, Eigen::aligned_allocator<VectorDIM>>;
    using Constraint = splx::Constraint<T>;

    PiecewiseCurveQPGenerator() : m_problem(0) {}

    void addPiece(std::shared_ptr<_QPOperations> opt_ptr) { //QPoperations
        m_operations.push_back(opt_ptr); // private member m_operations that is of vector type, thats why pushiing back operations pointer
        
        //Max_Param
        if(m_cumulativeMaxParameters.empty()) {
            m_cumulativeMaxParameters.push_back(opt_ptr->maxParameter());
        } else {
            m_cumulativeMaxParameters.push_back(
                m_cumulativeMaxParameters.back() + opt_ptr->maxParameter()
            );
        }

        // cummulative max params is filled using qpoperators...
        // max_param and numDecisionVariables are both included inside opt_ptr...
        // need to run a test case for this shit..


        //DecisionVar
        if(m_cumulativeDecisionVars.empty()) {
            m_cumulativeDecisionVars.push_back(opt_ptr->numDecisionVariables());
        } else {
            m_cumulativeDecisionVars.push_back(
                m_cumulativeDecisionVars.back() + opt_ptr->numDecisionVariables()
            );
        }

        m_problem = QPWrappers::Problem<T>(this->numDecisionVariables()); // referenced from qpoperators and appended to arg variable, stopped here...
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

    void addBezier(Index ncpts, T a) { // this was added in dyn_simulation solely.
        auto bezptr = std::make_shared<_BezierQPOperations>(ncpts, a);
        auto optptr = std::static_pointer_cast<_QPOperations>(bezptr); // returns a copy of the bezptr which is of type QPOperations..
        this->addPiece(optptr);
    }

    void setBezier(std::size_t idx, Index ncpts, T a) {
        auto bezptr = std::make_shared<_BezierQPOperations>(ncpts, a);
        auto optptr = std::static_pointer_cast<_QPOperations>(bezptr);
        this->setPiece(idx, optptr);
    }

    std::size_t numPieces() const {
        return m_cumulativeDecisionVars.size(); // 4 pieces of 8 ctrl pts each that are each of DIM D, every piece that is added, push back the number of ctrl_pts into decvar vector   
    }

    std::vector<T> pieceMaxParameters() const {
        std::vector<T> piece_max_params;
        for(const auto piece_opt : m_operations) {
            piece_max_params.push_back(piece_opt->maxParameter());
        }
        return piece_max_params;
    }

    T maxParameter() const {
        if(m_operations.empty()) {
            throw std::domain_error(
                absl::StrCat(
                    "piecewise curve has no pieces"
                )
            );
        }

        return m_cumulativeMaxParameters.back();
    }

    void setPieceMaxParameters(const std::vector<T>& new_max_params) {
        if(new_max_params.size() != m_operations.size()) {
            throw std::domain_error(
                absl::StrCat(
                    "new max parameters size does not match number of pieces",
                    ", number of pieces: ",
                    m_operations.size(),
                    ", given number of max parameters: ",
                    new_max_params.size()
                )
            );
        }

        for(std::size_t i = 0; i < m_operations.size(); i++) {
            m_operations[i]->maxParameter(new_max_params[i]);
        }

        this->fixCumulativeStructures(0);
        m_problem = QPWrappers::Problem<T>(this->numDecisionVariables());
    }

    Index numDecisionVariables() const {
        if(m_cumulativeDecisionVars.empty()) return 0; 

        return m_cumulativeDecisionVars.back(); // 24 cos its jus ctrl_pts * DIM, its back = 24, not size = 4
    } 
    
    void addIntegratedSquaredDerivativeCost(unsigned int k, T lambda) { // m_operations is of type QPoperations again ..
        Matrix Q(this->numDecisionVariables(), this->numDecisionVariables()); // 24
        Vector c(this->numDecisionVariables()); // 24
        Q.setZero();
        c.setZero();

        for(std::size_t i = 0; i < this->numPieces(); i++) { // 4
            Index dvar_start_idx = (i == 0 ? 0 : m_cumulativeDecisionVars[i-1]);
            Index dvar_count = m_operations[i]->numDecisionVariables();

            auto [Qs, cs]
                = m_operations[i]->integrantedSquaredDerivativeCost(k, lambda);
            
            Q.block(dvar_start_idx, dvar_start_idx, dvar_count, dvar_count) = Qs;
            c.block(dvar_start_idx, 0, dvar_count, 1) = cs;
        }

        m_problem.add_Q(Q);
        m_problem.add_c(c);
    }

    void addEvalCost(T u, unsigned int k, const VectorDIM& target, T lambda) {
        Matrix Q(this->numDecisionVariables(), this->numDecisionVariables());
        Vector c(this->numDecisionVariables());
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

    void addEvalConstraint(T u, unsigned int k, const VectorDIM& target,
                           bool soft_convertible = false,
                           T soft_convertible_weight = T(1)) {
        auto [idx, param, first_dvar_index, dvar_count] = pieceInfo(u);
        auto constraints = m_operations[idx]->evalConstraint(param, k, target,
                soft_convertible, soft_convertible_weight);

        this->addConstraints(constraints, first_dvar_index, dvar_count);
    }

    void addHyperplaneConstraintForPiece(std::size_t idx, const Hyperplane& hp,
                                         bool soft_convertible = false,
                                         T soft_convertible_weight = T(1)) {
        Index first_dvar_index = 
                (idx == 0 ? 0 : m_cumulativeDecisionVars[idx-1]);
        Index dvar_count = m_operations[idx]->numDecisionVariables();

        auto constraints = m_operations[idx]->hyperplaneConstraintAll(hp,
                                  soft_convertible, soft_convertible_weight);

        this->addConstraints(constraints, first_dvar_index, dvar_count);
    }

    void addHyperplaneConstraintAll(const Hyperplane& hp,
                                    bool soft_convertible = false,
                                    T soft_convertible_weight = T(1)) {
        for(std::size_t i = 0; i < this->numPieces(); i++) {
            this->addHyperplaneConstraintForPiece(i, hp,
                                  soft_convertible, soft_convertible_weight);
        }
    }

    void addHyperplaneConstraintAt(T u, const Hyperplane& hp,
                                   bool soft_convertible = false,
                                   T soft_convertible_weight = T(1)) {
        auto [idx, param, first_dvar_index, dvar_count] = pieceInfo(u);
        auto constraints = m_operations[idx]->hyperplaneConstraintAt(param, hp,
                                 soft_convertible, soft_convertible_weight);
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
    void addContinuityConstraint(std::size_t idx, unsigned int k,
                                 bool soft_convertible = false,
                                 T soft_convertible_weight = T(1)) {
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
                            d, 0, k
            );

            Row coeff(this->numDecisionVariables());
            coeff.setZero();
            coeff.block(0, first_piece_dvars_start, 1, first_piece_numdvars) 
                                = coeff1;
            coeff.block(0, second_piece_dvars_start, 1, second_piece_numdvars)
                                = -coeff2;

            m_problem.add_constraint(coeff, 0, 0, soft_convertible, soft_convertible_weight);
        }
    }

    Vector getDVarsForSegments(const StdVectorVectorDIM& segments) const {
        if(segments.size() != this->numPieces() + 1) {
            throw std::domain_error
            (
                absl::StrCat
                (
                    "number of segments is not equal to the number of pieces",
                    ", segment count: ",
                    segments.size() - 1,
                    ", piece count: ",
                    this->numPieces()
                )
            );
        }

        Vector res(this->numDecisionVariables());

        for(std::size_t i = 0; i < segments.size() - 1; i++) {
            Index piece_dvars_start
                = (i == 0 ? 0: m_cumulativeDecisionVars[i-1]);
            Index piece_numdvars = m_operations[i]->numDecisionVariables();
            res.block(piece_dvars_start, 0, piece_numdvars, 1)
                = m_operations[i]->getDVarsForSegment
            (
                        segments[i],
                        segments[i+1]
            );
        }

        return res;
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


    void resetProblem() {
        m_problem.reset();
    }

    void resetGenerator() {
        m_operations.clear();
        m_cumulativeMaxParameters.clear();
        m_cumulativeDecisionVars.clear();
        this->resetProblem();
    }

    const QPWrappers::Problem<T>& getProblem() const {
        return this->m_problem;
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

        return {
            idx,
            std::min(param, m_operations[idx]->maxParameter()),
            first_dvar_index,
            dvar_count
        };
    }

    void addConstraints(const std::vector<Constraint>& cons,
                        Index first_dvar_index, Index dvar_count) {

        Row coeff(this->numDecisionVariables());

        for(const auto& constraint: cons) {
            coeff.setZero();
            coeff.block(0, first_dvar_index, 1, dvar_count) 
                        = constraint.coeff;
            m_problem.add_constraint(coeff, 
                                     constraint.lb,
                                     constraint.ub,
                                     constraint.soft_convertible,
                                     constraint.soft_weight
            );   
        }
    }


};

}

#endif