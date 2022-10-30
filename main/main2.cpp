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
#include "pendulum_parametric.h"

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
#include "pendulum_parametric.h"
// ainda nao esta a funcionar :(

void Uma_Ode(); // Função para gerar uma só ode e guardar o output num ficheiro txt
int Is_Stable(double *, int);
std::vector<double> Moving_Average(std::vector<ODEpoint>, double, double);
void Stability_Map(int, double); // Gerar n ODE's, só quero testar umas coisas

int main()
{
    // Uma_Ode();
    Stability_Map(20, 0.1);

    return 0;
}

// Criterio de estabilidade do projeto
int Is_Stable(double *media, int n)
{
    int flag_estab = 1;
    for (int i = 0; i < n; i++)
    {
        if (media[i] > M_PI)
        {
            flag_estab = 0;
            break;
        }
    }
    return flag_estab;
};

std::vector<double> Moving_Average(std::vector<ODEpoint> vetor, double time_step, double tw)
{
    double sum_sq;
    std::vector<double> amp_mean;
    amp_mean.resize(vetor.size() - tw);
    int loop;

    for (int i = 0; i < vetor.size() - tw; i++)
    {
        amp_mean[i] = 0;
    }

    loop = 0;
    while (loop < tw)
    {
        sum_sq += vetor[loop].X()[0] * vetor[loop].X()[0];
        loop++;
    }
    double media = sqrt(sum_sq / tw);
    amp_mean[0] = media;

    for (int i = 1; i < vetor.size() - tw; i++)
    {
        sum_sq += vetor[i + tw].X()[0] * vetor[i + tw].X()[0];
        sum_sq -= vetor[i - 1].X()[0] * vetor[i - 1].X()[0];
        media = sqrt(sum_sq / tw);
        amp_mean[i] = media;
    }

    return amp_mean;
};

void Uma_Ode()
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
}

void Stability_Map(int numero, double step)
{
    int n;
    int count = 0;

    // Max and min delta values
    double delta_range[2] = {0.2, 0.35};
    double epsilon_range[2] = {0.075, 0.25};

    double step_delta = (delta_range[1] - delta_range[0]) / numero;
    double step_epsilon = (epsilon_range[1] - epsilon_range[0]) / numero;

    std::ofstream outdata;
    outdata.open("main/hist_data.txt"); // opens the file
    if (!outdata)
    { // file couldn't be openedcsv
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }

    double est = 1;

    std::cout << "A cozinhar esses gráficos... " << std::endl;

    for (double delta = delta_range[0]; delta <= delta_range[1]; delta = delta + step_delta)
    {
        // Adicionei loop unrolling aqui, daí andar de 2 em 2
        for (double epsilon = epsilon_range[0]; epsilon <= epsilon_range[1] + step_epsilon; epsilon = epsilon + 2 * step_epsilon)
        {
            pendulum_parametric pendulum(1, {0.1, 0});
            pendulum.SetFunction(0, [](ODEpoint p)
                                 { return p.X()[1]; });

            pendulum.SetFunction(1, [&](ODEpoint p)
                                 { return -((delta + epsilon * std::cos(p.T())) * std::sin((p.X())[0])); });

            std::vector<ODEpoint> resultado = pendulum.LeapFrogImprovedSolver(numero, step);
            std::vector<double> mov_avg = Moving_Average(resultado, step, 5 / step);

            n = resultado.size();

            double *media = &mov_avg[0];

            for (int i = 0; i < n; ++i)
            {
                outdata << resultado[i].X()[0];
                if (i != n - 1)
                    outdata << ";";
            }
            outdata << std::endl;
        }
        count++;
        if (count % 5 == 0)
            std::cout << "\u001b[32m" << 100 * (delta - delta_range[0]) / (delta_range[1] - delta_range[0]) << "% " << std::endl;
    }
    std::cout << " " << std::endl;

    outdata.close();

    // Correr o ficheiro python para fazer o dynamic mode decomposition
    /* std::cout << "\n\u001b[37mA correr o código python..." << std::endl;
    system("python3 main/Hist.py"); */
}