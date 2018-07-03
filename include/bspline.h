#ifndef SPLX_BSPLINE_H
#define SPLX_BSPLINE_H
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#include <vector>
#include "spline.h"
#include <utility>

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
     * Adds -alpha * (hp.normal (dot) pt - d) ^ 2 for each point from "from" to "to".
     * Effectively punishes the points as they get closer to the hyperplane
     *
     * order of variables is assumed to be
     * p0x, p1x, p2x, ..., pnx, p1y, ..., pny, ...
    */
    void extendQPHyperplanePenalty(QPMatrices& QP, unsigned int from, unsigned int to, const Hyperplane& hp, double alpha) const;

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
     * Add constraints that requires control points from index from to index to to be on the
     * negative side of hp.
    */
    void extendQPHyperplaneConstraint(QPMatrices& QP, unsigned int from, unsigned int to, const Hyperplane& hp) const;

    /**
     * Add constraint that requires curve to be on the negative side of hp
     * when from <= u <= to
     *
     * Effectively, enforces all points that effects the curve in [from, to] to be
     * in the negative side of the hp.
    */
    void extendQPHyperplaneConstraint(QPMatrices& QP, double from, double to, const Hyperplane& hp) const;

    /**
     * Add constraints that requires all decision variables to be less than ub and more than lb.
     *
     * Effectively, it constraints the curve to be inside a box.
    */
    void extendQPDecisionConstraint(QPMatrices&QP, double lb, double ub) const;
    /*
     * returns the inclusive index range of points that effects the curve
     * from u = from to u = to.
    */
    std::pair<unsigned int, unsigned int> affectingPoints(double from, double to) const;


    /*
     * Getter for k^th control point
    */
    const Vec& getCP(unsigned int k) const;

    /*
     * Interpolate from point 'from' to point 'to' with n points.
     * 'from' is the 1st point, 'to' is nth point.
     * Repeat 'to' m_degree+1 times to make the curve end at 'to'.
     * In the end there are n + m_degree points added.
     *
     * returns inclusive index range of added points in m_controlPoints array.
    */
    std::pair<unsigned int, unsigned int> interpolateEndAtTo(const Vec& from, const Vec& to, unsigned int n);


    /**
      * Generates knot vector from scratch
      *
      * @assumes control points are already set
    */
    void generateUniformKnotVector();

    /**
     * Clear control points array.
    */
    void clearControlPoints();

    /**
     * Load control points from QP.x
    */
    void loadControlPoints(const QPMatrices& QP);

    /**
      DBG FUNCTIONS
    */
    Vec eval_dbg(double u) const;
    void printKnotVector() const;
    void printControlPoints() const;
    void printKnotVectorNumbered() const;


    unsigned int m_degree; // degree of basis functions
    double m_a; // first p+1 knot values
    double m_b; // last p+1 knot values
    std::vector<Vec> m_controlPoints; // control points

  private:
    /**
      dimension that spline is defined in,
      defines the dimension of control points as well.
    */
    unsigned int m_dimension;
    std::vector<double> m_knotVector; // knot vector


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
