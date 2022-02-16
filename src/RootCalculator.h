#ifndef __RootCalculator__
#define __RootCalculator__

#include <cstdio>
#include <algorithm>
#include <stdexcept>
#include <iomanip> 
#include <iostream>
#include <vector>
#include <fstream>
#include <cstring>
#include <sstream>
#include <ostream>
#include <string>

#include "TCanvas.h" // include class TCanvas from ROOT library
#include "TRootCanvas.h"
#include "TH2F.h" // histogram 2D
#include "TApplication.h"
#include "TGraph.h"
#include "TF1.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TExec.h"
#include "TPolyMarker.h"

class RootCalculator
{
public:
    // constructors, destructor
    RootCalculator();
    RootCalculator(std::string expressao, double a,double b);
    ~RootCalculator() = default;  

    double firstDerivative(double x, double h, std::string option);
    std::vector<double> Secant_Method();

    void Draw(std::vector<double> zeros);

private:
    TFormula Funcao;
    std::string Expressao;
    double A,B;
};

#endif
