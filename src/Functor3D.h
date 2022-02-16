#ifndef __FUNCTOR3D__
#define __FUNCTOR3D__

#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

#include "TGraph2D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TAxis.h"
#include "TROOT.h"
#include "TH2F.h"
#include "TApplication.h"
#include "TRootCanvas.h"
#include <TROOT.h>
#include <TStyle.h>

class Functor3D
{
public:
    Functor3D(std::string s = "Functor3D") : name(s) { ; }
    ~Functor3D() = default;
    virtual double operator()(double x, double y);

    virtual void Draw(double xi, double xf, double yi, double yf, int num, std::string xtitle = "x", std::string ytitle = "y", std::string ztitle = "z");

protected:
    TCanvas *c;
    std::string name;
};

#endif
