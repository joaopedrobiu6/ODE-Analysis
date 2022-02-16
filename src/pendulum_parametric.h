#ifndef __PENDULUM_PARAMETRIC__
#define __PENDULUM_PARAMETRIC__

#include "ODEpoint.h"

#include "TGraph.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TAxis.h"
#include "TROOT.h"
#include "TH2F.h"
#include "TApplication.h"
#include "TRootCanvas.h"
#include <TROOT.h>
#include <TStyle.h>
#include "TMultiGraph.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TExec.h"
#include "TPolyMarker.h"
#include "TMath.h"
#include "TH2F.h"

#include <vector>
#include <iostream>
#include <sstream> // std::stringstream
#include <fstream> // std::fstream
#include <iomanip>
#include <limits>
#include <unistd.h>

class pendulum_parametric
{
public:
    // Constructor
    pendulum_parametric() = default;
    pendulum_parametric(double, const std::initializer_list<double> &); // pendulum_parametric(10, {80, 0})
    // length in meters, angles in degrees, velocity degrees/sec
    ~pendulum_parametric() = default;

    void SetFunction(int, std::function<double(ODEpoint)>);
    // solvers
    const std::vector<ODEpoint> &RungeKutta4(double Time = 0, double step = 1E-3);
    const std::vector<ODEpoint> &LeapFrogImprovedSolver(double Time = 0, double step = 1E-3); // método usado para confirmar valores do RK4

private:
    double L; // length (m)
    Xvar X0;  // initial conditions: angle (rad), angular velocity (rad/s)
    std::vector<ODEpoint> resultado;
    std::function<double(ODEpoint)> f[2]; // Recebe uma odepoint e retorna um double
};

#endif