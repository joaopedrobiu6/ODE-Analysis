#ifndef __DATAREADER__
#define __DATAREADER__

#include <iostream>
#include <vector>
#include <map>
#include <utility>
#include <cmath>
#include <sstream> // std::stringstream
#include <fstream> // std::fstream
#include <iomanip>

#include "TCanvas.h"
#include "TRootCanvas.h"
#include "TApplication.h"
#include "TGraph.h"
#include "TF1.h"
#include "TMultiGraph.h"

class DataReader
{
public:
    DataReader(std::string);
    DataReader(const DataReader &);
    ~DataReader() = default;

    const std::vector<std::vector<float>> &GetData();

    void dump();
    friend std::ostream &operator<<(std::ostream &, const DataReader &);

protected:
    std::vector<std::vector<float>> data;
};

#endif