#include <iostream>
#include <vector>
#include <functional>
#include <map>
#include <cmath>
#include <time.h>

#include "ODE_analysis.h"

int main()
{
    srand(time(NULL));
    double f = (double)rand() / RAND_MAX;

    double epsilon = (double)(rand() % 100) / 100;
    double delta = (double)(rand() % 100) / 100;
    std::cout << "(" << delta << ", " << epsilon << ")" << std::endl;

    double initial_position = M_PI / 2;
    ODE_analysis pendulum(2, {initial_position, 0});

    pendulum.SetFunction(0, [](ODEpoint p)
                         { return p.X()[1]; });

    pendulum.SetFunction(1, [&](ODEpoint p)
                         { return -((delta + epsilon * std::cos(p.T())) * std::sin((p.X())[0])); });

    std::vector<ODEpoint> result_euler = pendulum.EulerSolver(120, 0.0125);
    std::vector<ODEpoint> result_RK4 = pendulum.RungeKutta4(120, 1E-3);
    std::vector<ODEpoint> result_leapfrog = pendulum.LeapFrogImprovedSolver(120, 1E-3);

    pendulum.ODE_Draw("euler.png", "time", "#hat{#theta}", 50, result_euler);
    pendulum.ODE_Draw("RK4.png", "time", "#hat{#theta}", 50, result_RK4);
    pendulum.ODE_Draw("leapfrog.png", "time", "#hat{#theta}", 50, result_leapfrog);
    return 0;
}