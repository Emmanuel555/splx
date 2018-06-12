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
      * Evaluates spline at u
      *
      * @fails if u is not in [m_a, m_b]
    */
    Vec eval(double u);

    /**
      * Evaluates the k^{th} derivative of the spline at u \in [0,1]
      *
      * @fails if u is not in [m_a, m_b]
    */
    Vec eval(double u, unsigned int k);

    /**
      DBG FUNCTIONS
    */
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

  };

}

#endif
