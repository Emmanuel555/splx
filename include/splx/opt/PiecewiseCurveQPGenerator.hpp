#ifndef SPLX_PIECEWISECURVE_QP_GENERATOR_HPP
#define SPLX_PIECEWISECURVE_QP_GENERATOR_HPP

#include <qp_wrappers/problem.hpp>
#include <splx/opt/QPGenerator.hpp>
#include <splx/opt/BezierQPGenerator.hpp>
#include <splx/types.hpp>
#include <splx/curve/PiecewiseCurve.hpp>

namespace splx {

template<typename T, unsigned int DIM>
class PiecewiseCurveQPGenerator {
    public:
        using _QPGenerator = QPGenerator<T, DIM>;
        using _BezierQPGenerator = BezierQPGenerator<T, DIM>;
        using _Problem = QPWrappers::Problem<T>;
        using Row = splx::Row<T>;
        using VectorDIM = splx::VectorDIM<T, DIM>;
        using Hyperplane = splx::Hyperplane<T, DIM>;

        PiecewiseCurveQPGenerator() : m_problem(0), m_state(GeneratorState::CONFIGURE_PIECES) {}

        enum class GeneratorState {
            CONFIGURE_PIECES,
            FORMULATION,
            FORMULATION_DONE
        };

        void beginFormulation() {
            this->ensureStateIs(GeneratorState::CONFIGURE_PIECES);
            
            m_problem = _Problem(this->numControlPoints() * DIM);
            this->m_state = GeneratorState::FORMULATION;
        }

        void endFormulation() {
            this->ensureStateIs(GeneratorState::FORMULATION);

            Index vars_before = 0;
            Row constraint_coeff(m_problem.num_vars());
            for(unsigned int i = 0; i < m_subQPGenerators.size(); i++) {
                const _QPGenerator& subGenerator = *(m_subQPGenerators[i]);
                const _Problem& subProblem = subGenerator.getProblem();

                // merge cost
                m_problem.add_Q_block(vars_before, vars_before, subProblem.Q());
                m_problem.add_c_block(vars_before, subProblem.c());

                // merge var limits
                for(Index j = 0; j < subProblem.num_vars(); j++) {
                    m_problem.set_var_limits(vars_before + j, subProblem.lbx()(j), subProblem.ubx()(j));
                }

                // merge constraints
                for(Index j = 0; j < subProblem.num_constrainst(); j++) {
                    constraint_coeff.setZero();
                    constraint_coeff.block(0, vars_before, 1, subProblem.num_vars()) = subProblem.A().row(j);
                    m_problem.add_constraint(constraint_coeff, subProblem.lb()(j), subProblem.ub()(j));
                }

                vars_before += subGenerator.getProblem().num_vars();
            }

            this->m_state = GeneratorState::FORMULATION_DONE;
        }

        void resetFormulation() {
            this->ensureStateIs({GeneratorState::FORMULATION_DONE, GeneratorState::FORMULATION_DONE});

            for(auto subGenerator: m_subQPGenerators) {
                subGenerator->reset();
            }

            this->m_state = GeneratorState::CONFIGURE_PIECES;
        }


        /*
        * Add a bezier curve with ncpts control points and max parameter a
        */
        void addBezier(std::size_t ncpts, T a) {
            this->ensureStateIs(GeneratorState::CONFIGURE_PIECES);

            std::shared_ptr<_BezierQPGenerator> bezgen = std::make_shared<_BezierQPGenerator>(ncpts, a);
            m_subQPGenerators.push_back(std::static_pointer_cast<QPGenerator>(bezgen));
            m_cumulativeParameters.push_back(a);
        }

        void setBezier(std::size_t idx, std::size_t ncpts, T a) {
            this->ensureStateIs(GeneratorState::CONFIGURE_PIECES);
            this->pieceIndexCheck(idx);

            auto bezptr = std::make_shared<_BezierQPGenerator>(ncpts, a);
            auto genptr = std::static_pointer_cast<_QPGenerator>(bezptr);

            m_subQPGenerators[idx] = genptr;
            this->fixCumulativeParameters(idx);
        }

        void removePiece(std::size_t idx) {
            this->ensureStateIs(GeneratorState::CONFIGURE_PIECES);
            this->pieceIndexCheck(idx);

            m_subQPGenerators.erase(idx);
            m_cumulativeParameters.erase(idx);
            this->fixCumulativeParameters(idx);
        }

        void clearPieces() {
            this->ensureStateIs(GeneratorState::CONFIGURE_PIECES);
        }

        void addIntegratedSquaredDerivativeCost(unsigned int k, T lambda) {
            this->ensureStateIs(GeneratorState::FORMULATION);

            for(auto subGenerator : m_subQPGenerators) {
                subGenerator->addIntegratedSquaredDerivativeCost(k, lambda);
            }
        }

        void addIntegratedSquaredDerivativeCost(unsigned int k, const std::vector<T>& lambdas) {
            this->ensureStateIs(GeneratorState::FORMULATION);
            this->numPiecesCheck(lambdas.size());

            for(std::size_t i = 0; i < lambdas.size(); i++) {
                m_subQPGenerators[i]->addIntegratedSquaredDerivativeCost(k, lambdas[i]);
            }
        }

        void addEvalCost(T u, unsigned int k, const VectorDIM& target, T lambda) {
            this->ensureStateIs(GeneratorState::FORMULATION);
            auto [idx, ufixed] = this->pieceInfo(u);

            m_subQPGenerators[idx]->addEvalCost(ufixed, k, target, lambda);
        }

        void addEvalConstraint(T u, unsigned int k, const VectorDIM& target) {
            this->ensureStateIs(GeneratorState::FORMULATION);
            auto [idx, ufixed] = this->pieceInfo(u);

            m_subQPGenerators[idx]->addEvalConstraint(ufixed, k, target);
        }


        void addHyperplaneConstraintAll(const Hyperplane& hp) {
            this->ensureStateIs(GeneratorState::FORMULATION);
            for(auto subGenerator : m_subQPGenerators) {
                subGenerator->addHyperplaneConstraintAll(hp);
            }
        }


        void addHyperplaneConstraintAll(std::size_t idx, const Hyperplane& hp) {
            this->ensureStateIs(GeneratorState::FORMULATION);
            this->pieceIndexCheck(idx);

            m_subQPGenerators[idx]->addHyperplaneConstraintAll(hp);
        }


        void addHyperplaneConstraintAt(T u, const Hyperplane& hp) {
            this->ensureStateIs(GeneratorState::FORMULATION);
            auto [idx, ufixed] = this->pieceInfo(u);

            m_subQPGenerators[idx]->addHyperplaneConstraintAt(ufixed, hp);
        }


        void addControlPointLimits(const VectorDIM& lb, const VectorDIM& ub) {
            this->ensureStateIs(GeneratorState::FORMULATION);
            for(auto subGenerator : m_subQPGenerators) {
                subGenerator->addControlPointLimits(lb, ub);
            }
        }

        void addControlPointLimits(std::size_t idx, const VectorDIM& lb, const VectorDIM& ub) {
            this->ensureStateIs(GeneratorState::FORMULATION);
            this->pieceIndexCheck(idx);
            m_subQPGenerators[idx]->addControlPointLimits(lb, ub);
        }

        /*
        * continuity constraint between piece idx and piece idx+1
        */
        void addContinuityConstraint(std::size_t idx, unsigned int k) {
            this->ensureStateIs(GeneratorState::FORMULATION);
            this->pieceIndexCheck(idx);
            this->pieceIndexCheck(idx+1);

            Index vars_before = 0;
            for(std::size_t i = 0; i < idx; i++) {
                vars_before += m_subQPGenerators[i]->getProblem().num_vars();
            }


            Row coeff(this->numControlPoints() * DIM);

            for(unsigned int d = 0; d < DIM; d++) {
                Row coeff1 = m_subQPGenerators[idx]->evalBasisRow(
                    d, 
                    m_subQPGenerators[idx]->maxParameter(),
                    k
                );

                Row coeff2 = m_subQPGenerators[idx+1]->evalBasisRow(
                    d,
                    0,
                    k
                );

                coeff.setZero();
                coeff.block(0, vars_before, 1, coeff1.cols()) = coeff1;
                coeff.block(0, vars_before + coeff1.cols(), 1, coeff2.cols()) = -coeff2;

                m_problem.add_constraint(coeff, 0, 0);
            }
        }

        const _Problem& getProblem() const {
            this->ensureStateIs(GeneratorState::FORMULATION_DONE);
            return m_problem;
        }

        std::size_t numControlPoints() const {
            std::size_t ncpts = 0;
            for(const auto subGenerator : m_subQPGenerators) {
                ncpts += subGenerator->numControlPoints();
            }
            return ncpts;
        }

        T maxParameter() const {
            return m_cumulativeParameters.back();
        }

    private:
        std::vector<std::shared_ptr<_QPGenerator>> m_subQPGenerators;
        std::vector<T> m_cumulativeParameters;
        _Problem m_problem;
        GeneratorState m_state;

        std::pair<std::size_t, T> pieceInfo(T u) {
            this->parameterBoundCheck(u);
            std::size_t idx = std::lower_bound(
                    m_cumulativeParameters.begin(), 
                    m_cumulativeParameters.end(), 
                    u
                    ) - m_cumulativeParameters.begin();

            if(idx != 0)
                u -= m_cumulativeParameters[idx-1];

            return std::make_pair(idx, u);
        }

        void addNewCurveMaxParameter(T maxParam) {
            if(m_cumulativeParameters.empty()) {
                m_cumulativeParameters.push_back(maxParam);
            } else {
                m_cumulativeParameters.push_back(maxParam + m_cumulativeParameters.back());
            }
        }

        void fixCumulativeParameters(std::size_t idx) {
            for(; idx < m_cumulativeParameters.size(); idx++) {
                if(idx == 0) {
                    m_cumulativeParameters[idx] = m_subQPGenerators[idx]->maxParameter();
                } else {
                    m_cumulativeParameters[idx] = m_cumulativeParameters[idx-1] + m_subQPGenerators[idx]->maxParameter();
                }
            }
        }

        static std::string state_to_str(GeneratorState state) {
            switch(state) {
                case GeneratorState::CONFIGURE_PIECES:
                    return "CONFIGURE_PIECES";
                case GeneratorState::FORMULATION:
                    return "FORMULATION";
                case GeneratorState::FORMULATION_DONE:
                    return "FORMULATION_DONE";
                default:
                    throw std::domain_error("not handled state");
            }
        }

        void ensureStateIs(GeneratorState state) const {
            if(m_state != state) {
                std::string current_state = this->state_to_str(m_state), 
                            desired_state = this->state_to_str(state);
                throw std::logic_error(
                    std::string("current state is ")
                    + current_state
                    + std::string(" while desired state is ")
                    + desired_state
                );
            }
        }

        void ensureStateIs(const std::vector<GeneratorState>& states) const {
            if(std::find(states.begin(), states.end(), m_state) == states.end()) {
                throw std::logic_error(
                    std::string("state is not in the given set")
                );
            }
        }

        void pieceIndexCheck(std::size_t idx) const {
            if(idx >= m_cumulativeParameters.size()) {
                throw std::domain_error(
                    std::string("piece index out of range, given: ")
                    + std::to_string(idx)
                    + std::string(", num pieces: ")
                    + std::to_string(m_cumulativeParameters.size())
                );
            }
        }

        void numPiecesCheck(std::size_t num) const {
            if(num != m_cumulativeParameters.size()) {
                throw std::domain_error(
                    std::string("number of pieces does not match, given ")
                    + std::to_string(num)
                    + std::string(", have ")
                    + std::to_string(m_cumulativeParameters.size())
                );
            }
        }

        void emptyPiecesCheck() const { // checks if there is at least one piece.
            if(m_cumulativeParameters.empty()) {
                throw std::logic_error(
                    std::string("piecewise curve is empty.")
                );
            }
        }

        void parameterBoundCheck(T u) const { // checks if given parameter is valid
            this->emptyPiecesCheck();
            if(u < 0 || u > m_cumulativeParameters.back()) {
                throw std::domain_error(
                    std::string("given parameter is out of bounds. given u: ")
                    + std::to_string(u)
                    + std::string(", allowed range: [0, ")
                    + std::to_string(m_cumulativeParameters.back())
                    + std::string("]")
                );
        }
    }
};

}

#endif