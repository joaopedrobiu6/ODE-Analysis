#ifndef __Interpolator__
#define __Interpolator__

#include <iostream>
#include <vector>
#include <map>
#include <utility>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "DataPoints.h"
#include "EqSolver.h"
#include "FCmatrixAlgo.h"

#include "TCanvas.h" // include class TCanvas from ROOT library
#include "TRootCanvas.h"
#include "TApplication.h"
#include "TGraph.h"
#include "TF1.h"
#include "TMultiGraph.h"

class Interpolator : public DataPoints
{
public:
    // constructors, destructor
    Interpolator() = default;                                          // default constructor (nothing to be done?)
    Interpolator(int N = 0, double *x = nullptr, double *y = nullptr); // build DataPoints from C-arrays of x and yvalues
    Interpolator(const std::vector<std::pair<double, double>> &);
    Interpolator(const std::vector<double> &x, const std::vector<double> &y, double min, double max);
    Interpolator(const std::vector<std::pair<double, double>> &P1, double min, double max);

    ~Interpolator(); //=default;
    // interpolation methods
    // void Init();
    double InterpolateLagrange(double); // Lagrange interpolation
    void InterpolateNewtonCoefs();
    double InterpolateNewton(double);   // Newton interpolation
    double InterpolateSpline3(double);  // spline3
    // draw points and function
    void Draw(std::string s); // s="Lagrange", "Neville", "Spline3"
    // getters
    TF1 *GetFunction(std::string s); // s="Neville", "Spline3"

private:
    std::map<std::string, TF1 *> MI; // key="lagrange", "newton", "spline3"
    std::vector<double> Coefs;
};

#endif