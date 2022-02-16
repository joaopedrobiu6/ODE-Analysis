#ifndef __MATERIAL__
#define __MATERIAL__

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

class Material {
    public:
    Material(std::string fname="", Double_t fdens=0): name(fname), density(fdens) {;}
    std::string GetName() {return name;}
    Double_t GetDensity() {return density;}
    virtual void Print();

    protected:
    std::string name;
    Double_t density;
};

#endif