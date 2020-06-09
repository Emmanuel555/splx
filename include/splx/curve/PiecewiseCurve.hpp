#ifndef SPLX_PIECEWISECURVE_HPP
#define SPLX_PIECEWISECURVE_HPP
#include <memory>
#include <splx/curve/ParametricCurve.hpp>
#include <splx/curve/Bezier.hpp>

namespace splx {

template<typename T, unsigned int DIM>
class PiecewiseCurve {
public:
    using _ParametricCurve = splx::ParametricCurve<T, DIM>;
    using _Bezier = splx::Bezier<T, DIM>;
    using CurveType = typename _ParametricCurve::CurveType;
    using VectorDIM = typename _ParametricCurve::VectorDIM;

    PiecewiseCurve() {

    }

    ~PiecewiseCurve() {

    }

    /*
    * Adds a bezier piece to the curve by creating a copy.
    */
    void addPiece(const _Bezier& bez) {
        std::shared_ptr<_Bezier> bezptr = std::make_shared<_Bezier>(bez);
        std::shared_ptr<_ParametricCurve> paramptr = std::static_pointer_cast<_ParametricCurve>(bezptr);
        m_pieces.push_back(paramptr);

        this->addNewCurveMaxParameter(bez.maxParameter());
    }

    /*
    * Adds a bezier piece to the curve directly. The object must not be destroyed outside!
    */
    void addPiece(std::shared_ptr<_Bezier> bezptr) {
        std::shared_ptr<_ParametricCurve> paramptr =  std::static_pointer_cast<_ParametricCurve>(bezptr);
        m_pieces.push_back(paramptr);

        this->addNewCurveMaxParameter(bezptr->maxParameter());
    }

    // add a piece directly
    void addPiece(std::shared_ptr<_ParametricCurve> pieceptr) {
        m_pieces.push_back(pieceptr);
        this->addNewCurveMaxParameter(pieceptr->maxParameter());
    }

    /*
    * Sets a specific piece to the given bezier. Given bezier is copied
    */
    void setPiece(std::size_t idx, const _Bezier& bez) {
        this->pieceIndexCheck(idx);

        std::shared_ptr<_Bezier> bezptr = std::make_shared<_Bezier>(bez);
        std::shared_ptr<_ParametricCurve> paramptr = std::static_pointer_cast<_ParametricCurve>(bezptr);
        m_pieces[idx] = paramptr;

        this->fixCumulativeParameters(idx);
    }

    /*
    * Sets a specific piece to the given bezier. Given bezier must not be destroyed outside!
    */
    void setPiece(std::size_t idx, std::shared_ptr<_Bezier> bezptr) {
        this->pieceIndexCheck(idx);

        auto paramptr = std::static_pointer_cast<_ParametricCurve>(bezptr);
        m_pieces[idx] = paramptr;

        this->fixCumulativeParameters(idx);
    }

    void setPiece(std::size_t idx, std::shared_ptr<_ParametricCurve> pieceptr) {
        this->pieceIndexCheck(idx);
        m_pieces[idx] = pieceptr;
        this->fixCumulativeParameters(idx);
    }

    // get number of pieces
    std::size_t numPieces() const {
        return this->m_cumulativeParameters.size();
    }


    // get piece with given index
    const _Bezier& operator[](std::size_t idx) const {
        this->pieceIndexCheck(idx);

        if(m_pieces[idx]->type != CurveType::BEZIER) {
            throw std::domain_error(
                std::string("piece with index ")
                + std::to_string(idx)
                + std::string(" is not a bezier")
            );
        }

        auto bezptr = std::static_pointer_cast<_Bezier>(m_pieces[idx]);
        return *bezptr;
    }
    const _Bezier& getPiece(std::size_t idx) const {
        return this->operator[](idx);
    }

    // get the type of the piece with the given index
    CurveType type(std::size_t idx) const {
        this->pieceIndexCheck(idx);
        return m_pieces[idx]->type;
    }


    // evaluate the kth derivative at piecewise curve at parameter u
    VectorDIM eval(T u, unsigned int k) const {
        this->parameterBoundCheck(u);
        this->emptyPiecesCheck();

        auto idx = std::lower_bound(m_cumulativeParameters.begin(), m_cumulativeParameters.end(), u)
                   - m_cumulativeParameters.begin();


        if(idx != 0)
            u -= m_cumulativeParameters[idx-1];

        return m_pieces[idx]->eval(
                std::min(u, m_pieces[idx]->maxParameter()),
                k
        );
    }

    T maxParameter() const {
        this->emptyPiecesCheck();
        return m_cumulativeParameters.back();
    }

private:
    std::vector<std::shared_ptr<_ParametricCurve>> m_pieces;
    std::vector<T> m_cumulativeParameters;

    void addNewCurveMaxParameter(T maxParam) {
        if(m_cumulativeParameters.empty()) {
            m_cumulativeParameters.push_back(maxParam);
        } else {
            m_cumulativeParameters.push_back(m_cumulativeParameters.back() + maxParam);
        }
    }

    void fixCumulativeParameters(std::size_t idx) {
        for(; idx < m_cumulativeParameters.size(); idx++) {
            if(idx == 0) {
                m_cumulativeParameters[idx] = m_pieces[idx]->maxParameter();
            } else {
                m_cumulativeParameters[idx] = m_cumulativeParameters[idx-1] + m_pieces[idx]->maxParameter();
            }
        }
    }

    void pieceIndexCheck(std::size_t idx) const { // checks if piece index is valid
        if(idx >= m_cumulativeParameters.size() || idx < 0) {
            throw std::domain_error(
                std::string("piece index used is ")
                + std::to_string(idx)
                + std::string(" while piece count is ")
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

} // namespace splx

#endif