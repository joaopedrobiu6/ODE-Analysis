#include <iostream>
#include <vector>
#include <functional>
#include <map>
#include <cmath>
#include <time.h>
#include <fstream>
#include <cstdlib>

#include "ODE_analysis.h"
#include "tools.h"

// ainda nao esta a funcionar :(

int main()
{
    double initial_position = 0.1;
    ODE_analysis pendulum(2, {initial_position, 0});

    double delta = 0.1;
    double epsilon = 0.1;

    std::vector<ODEpoint> result_leapfrog;
    result_leapfrog.resize(150);

    int n = (int)result_leapfrog.size();

    std::ofstream outdata;
    outdata.open("data.csv"); // opens the file
    if (!outdata)
    { // file couldn't be openedcsv
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }

    for (int i = 0; i < 10; i++)
    {
        for (int j = 0; j < 10; j++)
        {
            pendulum.SetFunction(0, [](ODEpoint p)
                                 { return p.X()[1]; });

            pendulum.SetFunction(1, [&](ODEpoint p)
                                 { return -((delta + epsilon * std::cos(p.T())) * std::sin((p.X())[0])); });

            result_leapfrog = pendulum.LeapFrogImprovedSolver(15, 1e-1);

            for (int i = 0; i < n; ++i)
            {
                outdata << result_leapfrog[i].X()[0];
                if (i != n - 1)
                    outdata << ";";
            }
            outdata << std::endl;
            epsilon = epsilon + 0.1;
        }
        delta = delta + 0.1;
    }

    outdata.close();

    // Correr o ficheiro python para fazer o dynamic mode decomposition
    std::cout << "\nA correr o código python..." << std::endl;
    system("python3 main/ODEdmd.py");

    return 0;
}
