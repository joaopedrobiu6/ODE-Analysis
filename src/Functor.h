#ifndef __FUNCTOR__
#define __FUNCTOR__
 
#include <string>
#include <iostream>
#include <vector>
#include <cmath> 
#include <iomanip>

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


class Functor
{
public:
    Functor(std::string s = "Functor") : name(s) {;}
    ~Functor() = default;
    virtual double operator()(double x);
    // args:
    // xi, xf ........... xmin and xmax limits for function display
    // num .............. number of sampling points to be used o TGraph
    // xtitle, ytitle ... axis titles
    virtual void Draw(double xi, double xf, int num, std::string xtitle = "x", std::string ytitle = "y", std::string option = "");

protected:
    TCanvas *c;
    std::string name;
};

#endif