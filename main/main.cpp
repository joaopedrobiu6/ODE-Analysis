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

    // std::vector<ODEpoint> result_euler = pendulum.EulerSolver(120, 0.0125);
    //std::vector<ODEpoint> result_RK4 = pendulum.RungeKutta4(30, 1E-1);
    std::vector<ODEpoint> result_leapfrog = pendulum.LeapFrogImprovedSolver(30, 1E-1);

    // pendulum.ODE_Draw("images/euler.png", "time", "#hat{#theta}", 50, result_euler);
    // pendulum.ODE_Draw("images/RK42.png", "time", "#hat{#theta}", 50, result_RK4);
    // pendulum.ODE_Draw("images/leapfrog.png", "time", "#hat{#theta}", 50, result_leapfrog);
    // pendulum.Draw_Stability("images/stability.pdf", 150, 0.0125);

    std::ofstream outdata;
    outdata.open("data1.txt"); // opens the file
    if (!outdata)
    { // file couldn't be openedcsv
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }

    int n = (int)result_leapfrog.size();
    for (int i = 0; i < n; ++i)
    {
        outdata << result_leapfrog[i].T() << ";" << result_leapfrog[i].X()[0] << std::endl;
    }
    outdata.close();

    //Correr o ficheiro python para fazer o dynamic mode decomposition
    std::cout << "\nA correr o código python..." << std::endl;
    system("python3 main/ODEdmd.py");

    return 0;
}
