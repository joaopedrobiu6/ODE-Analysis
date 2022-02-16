#ifndef __FOURIERDT__
#define __FOURIERDT__

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

class FourierDT
{
public:
    FourierDT(std::vector<ODEpoint> data, double step_) : resultado(data), step(step_) { ; }
    ~FourierDT();

    std::vector<std::pair<double, double>> Fourier();
    void Fourier_Draw_Harmonics(const char *s, std::vector<std::pair<double, double>>);
    void Fourier_Draw_Max(const char *s, std::vector<std::pair<double, double>>, int qual_deles);
    void Comparacao_analitica(const char *filename, const char *xaxis, const char *yaxis, int color, std::vector<std::pair<double, double>> Freq);

private:
    double step;
    std::vector<ODEpoint> resultado;
};

#endif