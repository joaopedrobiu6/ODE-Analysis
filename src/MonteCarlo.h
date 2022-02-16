#ifndef __MONTECARLO__
#define __MONTECARLO__

#include "TAxis.h"
#include "TROOT.h"
#include "TF1.h"
#include "TH2F.h"
#include "TApplication.h"
#include "TRootCanvas.h"
#include <TROOT.h>
#include <TStyle.h>
#include "TH1F.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TDirectory.h"
#include "TFrame.h"
#include "TGraphErrors.h"
#include "TRandom3.h"

#include <vector>
#include <iostream>
#include <sstream> // std::stringstream
#include <fstream> // std::fstream
#include <iomanip>
#include <limits>
#include "Functor.h"
#include "MyFunction.h"

class MonteCarlo
{
public:
    MonteCarlo() = default;
    MonteCarlo(int N_, Functor &Func) : N(N_), F(Func) { ; }
    ~MonteCarlo() = default;

    // Getter Number of Samples
    int GetSampleNumber() { return N; };

    // MonteCarlo Integration
    void MCIntegration(double xi, double xf, double &Integral, double &Error);

private:
    int N;
    Functor &F;
};

#endif