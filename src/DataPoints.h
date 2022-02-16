#ifndef __DataPoints__
#define __DataPoints__

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

class DataPoints
{
public:
    // constructors, destructor
    DataPoints() = default;              // default constructor (nothing to be done?)
    DataPoints(int, double *, double *); // build DataPoints from C-arrays of x and y values
    DataPoints(const std::vector<std::pair<double, double>> &);
    DataPoints(const std::vector<double> &, const std::vector<double> &);
    ~DataPoints() = default;

    // getters
    const std::vector<std::pair<double, double>> &GetPoints();
    TGraph& GetGraphPointer();
    void GetGraph(TGraph &);

    // draw points using ROOT object TGraph
    virtual void Draw();

    // friend functions (optional)
    friend std::ostream &operator<<(std::ostream &, const DataPoints &);

protected:
    std::vector<std::pair<double, double>> P;
    TGraph *graph;
};

#endif