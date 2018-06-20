#ifndef SPLX_BSPLINE_H
#define SPLX_BSPLINE_H
#include <eigen3/Eigen/Dense>
#include <vector>
#include "spline.h"

namespace splx {

  class BSpline: public Spline {
  public:

    /**
     * Construct a b-spline in dim dimensions where basis functions are of
     * degree deg and curve has parameter u in [A, B]
     *
     * @fails if A > B
    */
    BSpline(unsigned int deg, unsigned int dim, double A, double B);

    /**
     * Construct a b-spline in dim dimensions where basis functions are of
     * degree deg and curve has parameter u \in [A, B] with given initial
     * control points
     *
     * @fails if A >= B
     * @fails if the size of m_controlPoints is less than m_degree+1
    */
    BSpline(unsigned int deg, unsigned int dim, double A,
            double B, const std::vector<Vec>& cpts);


    /**
      * Evaluates the k^{th} derivative of the spline at u \in [0,1]
      *
      * @fails if u is not in [m_a, m_b]
    */
    Vec eval(double u, unsigned int k) const;


    /**
     * In all functions that extends these matrices, order of variables is assumed to be
     * p0x, p1x, p2x, ..., pnx, p0y, p1y, ..., pny, ...
     *
    */
    QPMatrices getQPMatrices() const;


    /**
     * Add integral from m_a to m_b of square of norm of k^th derivative of the spline
     * to the hessian matrix H with scalar lambda
     *
     * order of variables is assumed to be
     * p0x, p1x, p2x, ..., pnx, p0y, p1y, ..., pny, ...
    */
    void extendQPIntegratedSquaredDerivative(QPMatrices& QP, unsigned int k, double lambda) const;


    /**
     * Add the cost theta * ||f(u) - pos||^2 to H and g.
     *
     * order of variables is assumed to be
     * p0x, p1x, p2x, ..., pnx, p0y, p1y, ..., pny, ...
    */

    void extendQPPositionAt(QPMatrices& QP, double u, const Vec& pos, double theta) const;

    /**
     * Add constraint that requires the k^th derivative of spline at u=0 to be target.
    */
    void extendQPBeginningConstraint(QPMatrices& QP, unsigned int k, const Vec& target) const;
    /**
      DBG FUNCTIONS
    */
    Vec eval_dbg(double u) const;
    void printKnotVector() const;
    void printControlPoints() const;
    void printKnotVectorNumbered() const;

  private:
    unsigned int m_degree; // degree of basis functions
    /**
      dimension that spline is defined in,
      defines the dimension of control points as well.
    */
    unsigned int m_dimension;
    double m_a; // first p+1 knot values
    double m_b; // last p+1 knot values
    std::vector<double> m_knotVector; // knot vector
    std::vector<Vec> m_controlPoints; // control points
    /**
      * Generates knot vector from scratch
      *
      * @assumes control points are already set
    */
    void generateUniformKnotVector();


    /**
      * Finds the index i of m_knotVector where u falls into [u_i, u_{i+1})
      *
      * @fails if u not in [0,1]
    */
    unsigned int findSpan(double u) const;

    /**
     * Evaluate k^th derivative of basis functions [N_{from,deg}, ..., N_{to,deg}] at u.
    */
    std::vector<double> evalBasisFuncs(double u, unsigned int deg, unsigned int k, unsigned int from, unsigned int to) const;

    /**
     * Get coefficient matrix of basis functions [N_{from,p}(u) ... N_{to,p}(u)] in interval
     * [m_knotVector[i], m_knotVector[i+1]) where first row is the coefficients of N_{from, p}
     * and last row is the coefficients of N_{to, p}
     *
     * a0 + a1u + a2u^2 + ... + apu^p
    */
    Matrix getBasisCoefficientMatrix(unsigned int from, unsigned int to, unsigned int p, unsigned int i) const;
  };

}

#endif
