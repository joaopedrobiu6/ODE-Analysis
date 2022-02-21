#ifndef __TOOLS__
#define __TOOLS__

#include <iostream>
#include <vector>
#include <functional>
#include <map>
#include <cmath>
#include <time.h>
#include <fstream>
#include <cstdlib>

#include "ODE_analysis.h"

void WriteData(const char *filename, std::vector<ODEpoint> P)
{
    std::ofstream outdata;
    outdata.open(filename); // opens the file
    if (!outdata)
    { // file couldn't be openedcsv
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }

    int n = (int)P.size();
    for (int i = 0; i < n; ++i)
    {
        outdata << P[i].T() << "," << P[i].X()[0] <<  std::endl;
    }
    outdata.close();
};

#endif
