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
    double delta = 0.7;
    double epsilon = 0.2;
    std::cout << "(" << delta << ", " << epsilon << ")" << std::endl;

    double initial_position = 0.1;
    ODE_analysis pendulum(2, {initial_position, 0});

    pendulum.SetFunction(0, [](ODEpoint p)
                         { return p.X()[1]; });

    pendulum.SetFunction(1, [&](ODEpoint p)
                         { return -((delta + epsilon * std::cos(p.T())) * std::sin((p.X())[0])); });

    std::vector<ODEpoint> result_RK4 = pendulum.RungeKutta4(30, 1E-1);

    WriteData("data.txt", result_RK4);

    // Correr o ficheiro python para fazer o dynamic mode decomposition
    std::cout << "\nA correr o código python..." << std::endl;
    system("python3 main/ODEdmd.py");

    return 0;
}
