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

int main()
{
    srand(time(NULL));
    double delta = (double)(rand() % 100) / 10;
    double epsilon = (double)(rand() % 100) / 10;
    std::cout << "(" << delta << ", " << epsilon << ")" << std::endl;

    double initial_position = 0.1;
    ODE_analysis pendulum(2, {initial_position, 0});

    pendulum.SetFunction(0, [](ODEpoint p)
                         { return p.X()[1]; });

    pendulum.SetFunction(1, [&](ODEpoint p)
                         { return -((delta + epsilon * std::cos(p.T())) * std::sin((p.X())[0])); });

    std::vector<ODEpoint> result_leapfrog = pendulum.LeapFrogImprovedSolver(120, 1e-1);

    WriteData("data1.txt", result_leapfrog);

    // Correr o ficheiro python para fazer o dynamic mode decomposition
    std::cout << "\nA correr o código python..." << std::endl;
    system("python3 main/ODEdmd.py");

    return 0;
}
