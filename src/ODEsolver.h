#ifndef __ODEsolver__
#define __ODEsolver__

#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <map>

#include "ODEpoint.h"

#include "TGraph.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TAxis.h"
#include "TROOT.h"
#include "TH2F.h"
#include "TApplication.h"
#include "TRootCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMultiGraph.h"

class ODEsolver
{
public:
    ODEsolver(const std::initializer_list<double> &X);
    ODEsolver() = default;
    ODEsolver(const Xvar &);
    ~ODEsolver() = default;

    void SetFunction(int, std::function<double(ODEpoint)>);
    void GetDim();

    const std::vector<ODEpoint> &EulerSolver(double Time = 0, double step = 1E-3);
    const std::vector<ODEpoint> &RungeKutta4(double Time = 0, double step = 1E-3);
    const std::vector<ODEpoint> &LeapFrogImprovedSolver(double Time = 0, double step = 1E-3);

private:
    Xvar X0; // initial conditions
    // solutions
    std::map<std::string, std::vector<ODEpoint>> MS; // key: "verlet", "euler","trapezoidal"
    // functions associated to dependent variables 1st order ODE's
    std::function<double(ODEpoint)> *f; // Recebe uma odepoint e retorna um double
};

#endif