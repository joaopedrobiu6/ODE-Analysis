#include "ODEsolver.h"

ODEsolver::ODEsolver(const std::initializer_list<double> &X)
{
    Xvar var(X);
    X0 = var;
};

ODEsolver::ODEsolver(const Xvar &X_)
{
    X0 = X_;
    int dim = X0.X().size();
    f = new std::function<double(ODEpoint)>[dim]();
};

void ODEsolver::SetFunction(int index, std::function<double(ODEpoint)> f_)
{
    f[index] = f_;
};

void ODEsolver::GetDim() { std::cout << X0.X().size() << std::endl; };

const std::vector<ODEpoint> &ODEsolver::EulerSolver(double Time, double step)
{
    // Time = tempo final
    // n  = quantos pontos vamos retirar entre tempo0 e Time
    /* basicamente em cada n (instante tempo0 < t < Time) vamos calcular um os Xvar (variaveis dependentes)
    (os instantes sao ODEpoints - a cada um deles esta associado um Xvar) */

    std::vector<ODEpoint> info1;
    int n = Time / step;

    // resizes
    info1.resize(n);        // precisamos que haja tantos ODEpoints quanto instantes de tempo!
    int n1 = X0.X().size(); // n1 é o tamanho de x de X0 - número de variáveis dependentes (dimensão dos Xvar!!!)

    for (int i = 0; i < n; i++)
        info1[i].X().resize(n1); // redimensionar o Xvar de cada um de n ODEpoints

    // condições iniciais - o primeiro ODEpoint que corresponde às condições do construtor!!!!
    info1[0].T() = 0;

    for (int i = 0; i < n1; i++)
        info1[0].X()[i] = X0.X()[i]; // colocar os dados do construtor no primeiro ODEpoint!!!

    // calculo das soluções

    for (int i = 0; i < n - 1; i++)
    {
        info1[i + 1].T() = step + info1[i].T();
        for (int j = 0; j < n1; j++)
        {
            info1[i + 1].X()[j] = info1[i].X()[j] + step * f[j](info1[i]);
        }
    }

    MS.insert({"euler", info1});
    return MS["euler"];
}

const std::vector<ODEpoint> &ODEsolver::RungeKutta4(double Time, double step)
{
    std::vector<ODEpoint> resultado;
    std::vector<double> K1, K2, K3, K4;

    // resizes
    int n = Time / step;
    int n1 = X0.X().size();

    resultado.resize(n);
    K1.resize(n1);
    K2.resize(n1);
    K3.resize(n1);
    K4.resize(n1);

    for (int i = 0; i < n; i++)
    {
        resultado[i].X().resize(n1);
    }

    // condições iniciais no tempo
    resultado[0].T() = 0;

    for (int i = 0; i < n1; i++)
    {
        resultado[0].X()[i] = X0.X()[i];
    }

    ODEpoint p0;
    Xvar xvar_temp;
    double tempo;
    p0.X().resize(n1);
    xvar_temp.X().resize(n1);

    for (int i = 0; i < n - 1; i++)
    {
        // Calcular o tempo em i+1
        resultado[i + 1].T() = resultado[i].T() + step;

        // Cálculo dos K1
        K1[0] = f[0](resultado[i]);
        K1[1] = f[1](resultado[i]);

        xvar_temp.X()[0] = resultado[i].X()[0] + (0.5 * step * K1[0]);
        xvar_temp.X()[1] = resultado[i].X()[1] + (0.5 * step * K1[1]);
        tempo = resultado[i].T() + (step * 0.5);

        p0.SetODEpoint(tempo, xvar_temp);

        // Cálculo dos K2
        K2[0] = f[0](p0);
        K2[1] = f[1](p0);

        xvar_temp.X()[0] = resultado[i].X()[0] + (0.5 * step * K2[0]);
        xvar_temp.X()[1] = resultado[i].X()[1] + (0.5 * step * K2[1]);

        p0.SetODEpoint(tempo, xvar_temp);

        // Cálculo dos K3
        K3[0] = f[0](p0);
        K3[1] = f[1](p0);

        // Cálculo dos K4
        xvar_temp.X()[0] = resultado[i].X()[0] + (step * K3[0]);
        xvar_temp.X()[1] = resultado[i].X()[1] + (step * K3[1]);
        tempo = resultado[i].T() + (step);

        p0.SetODEpoint(tempo, xvar_temp);

        K4[0] = f[0](p0);
        K4[1] = f[1](p0);

        // Cálculo das soluções
        resultado[i + 1].X()[0] = resultado[i].X()[0] + step * (K1[0] + 2 * K2[0] + 2 * K3[0] + K4[0]) / 6;
        resultado[i + 1].X()[1] = resultado[i].X()[1] + step * (K1[1] + 2 * K2[1] + 2 * K3[1] + K4[1]) / 6;

        K1.clear();
        K2.clear();
        K3.clear();
        K4.clear();
    }
    MS.insert({"rungekutta4", resultado});
    return MS["rungekutta4"];
};

const std::vector<ODEpoint> &ODEsolver::LeapFrogImprovedSolver(double Time, double step)
{
    std::vector<ODEpoint> info4;
    int n = Time / step;

    // resizes
    info4.resize(n);
    int n1 = X0.X().size(); // n1 é o tamanho de x de X0 - número de variáveis dependentes

    for (int i = 0; i < n; i++)
        info4[i].X().resize(n1);

    // condições iniciais
    info4[0].T() = 0;

    for (int i = 0; i < n1; i++)
        info4[0].X()[i] = X0.X()[i];

    // primeira iteração
    info4[1].T() = step + info4[0].T();
    for (int k = 0; k < n1; k++)
        info4[1].X()[k] = info4[0].X()[k] + step * f[k](info4[0]);

    for (int i = 1; i < n - 1; i++)
    {
        info4[i + 1].T() = step + info4[i].T();

        for (int k = 0; k < n1; k++)
        {
            info4[i + 1].X()[k] = info4[i].X()[k] + step * f[k](info4[i]) + 0.5 * step * step * f[k + 1](info4[i]);
            info4[i + 1].X()[k + 1] = info4[i].X()[k + 1] + 0.5 * step * (f[k + 1](info4[i]) + f[k + 1](info4[i + 1]));
        }
    }

    MS.insert({"leapfrog", info4});
    return MS["leapfrog"];
};