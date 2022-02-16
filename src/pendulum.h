#ifndef __PENDULUM__
#define __PENDULUM__

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

class pendulum
{
public:
    // Constructor
    pendulum() = default;
    pendulum(double, const Xvar &);
    pendulum(double, const std::initializer_list<double> &); // pendulum(10, {80, 0})
    // length in meters, angles in degrees, velocity degrees/sec
    pendulum(double length, double theta_0 = 80, double theta_vel_0 = 0);
    ~pendulum() = default;

    // solvers
    const std::vector<ODEpoint> &EulerSolver(double Time = 0, double step = 1E-3);
    const std::vector<ODEpoint> &EulerCromerSolver(double Time = 0, double step = 1E-3);
    const std::vector<ODEpoint> &VerletSolver(double Time = 0, double step = 1E-3);
    const std::vector<ODEpoint> &TrapezoidalSolver(double Time = 0, double step = 1E-3);
    const std::vector<ODEpoint> &LeapFrogImprovedSolver(double Time = 0, double step = 1E-3);
    const std::vector<ODEpoint> &RungeKutta2(double Time = 0, double step = 1E-3);
    const std::vector<ODEpoint> &RungeKutta4(double Time = 0, double step = 1E-3);

    // Draw
    void Draw(std::string s = "verlet", double ti = 0., double tf = 10.);

private:
    double L; // length (m)
    Xvar X0;  // initial conditions: angle (rad), angular velocity (rad/s)

    // solutions
    std::map<std::string, std::vector<ODEpoint>> MS; // key: "verlet", "euler","trapezoidal"

    // functions associated to dependent variables 1st order ODE's
    std::function<double(ODEpoint)> f[2]; // Recebe uma odepoint e retorna um double
};

#endif