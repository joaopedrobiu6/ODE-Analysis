#ifndef __RWALK1D__
#define __RWALK1D__

#include <iostream>
#include <algorithm>
#include <vector>
#include <iomanip>
#include <time.h>
#include <map>

#include "TCanvas.h" // include class TCanvas from ROOT library
#include "TRootCanvas.h"
#include "TH2F.h" // histogram 2D
#include "TApplication.h"
#include "TGraph.h"
#include "TF1.h"
#include "TROOT.h"
#include "TMultiGraph.h"

class Rwalk1D
{
public:
    Rwalk1D(int N_ = 1, double x_ = 0.,         // N=nb of particles, x=x(0)
            double pL_ = 0.5, double pR_ = 0.5, // probabilities Left Right
            double dt_ = 1, double dx_ = 1      // time and space steps
    );
    ~Rwalk1D();

    // particle simulation
    void Run(int nsteps); // number of steps

    // getters
    const std::vector<double>& GetTrajectory(int n = 1); // particle number
    double GetTimeStep();
    double GetSpaceStep();

private:
    double x0;                             // init coo
    int N;                                 // number of particles
    double pL, pR;                         // probabilities (left, same, righ)
    double dt, dx;                         // steps (time, space)
    std::map<int, std::vector<double>> mT; // trajectories
};

#endif