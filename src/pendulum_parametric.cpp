#include "pendulum_parametric.h"

pendulum_parametric::pendulum_parametric(double l, const std::initializer_list<double> &X)
{
    Xvar var(X);
    X0 = var;
    L = l;
};

void pendulum_parametric::SetFunction(int indice, std::function<double(ODEpoint)> f_)
{
    if (0 <= indice && indice <= 1)
        f[indice] = f_;
    else
        std::cout << "\nError, there can only be two equations (index 0 and 1)\n"
                  << std::endl;
};

const std::vector<ODEpoint> &pendulum_parametric::RungeKutta4(double Time, double step)
{
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

    return resultado;
};

const std::vector<ODEpoint> &pendulum_parametric::LeapFrogImprovedSolver(double Time, double step)
{
    // std::vector<ODEpoint> resultado;
    int n = Time / step;

    // resizes
    resultado.resize(n);
    int n1 = X0.X().size(); // n1 é o tamanho de x de X0 - número de variáveis dependentes

    for (int i = 0; i < n; i++)
        resultado[i].X().resize(n1);

    // condições iniciais
    resultado[0].T() = 0;

    for (int i = 0; i < n1; i++)
        resultado[0].X()[i] = X0.X()[i];

    // primeira iteração
    resultado[1].T() = step + resultado[0].T();
    resultado[1].X()[0] = resultado[0].X()[0] + step * f[0](resultado[0]);
    resultado[1].X()[1] = resultado[0].X()[1] + step * f[1](resultado[0]);

    for (int i = 1; i < n - 1; i++)
    {
        resultado[i + 1].T() = step + resultado[i].T();
        resultado[i + 1].X()[0] = resultado[i].X()[0] + step * f[0](resultado[i]) + 0.5 * step * step * f[1](resultado[i]);
        resultado[i + 1].X()[1] = resultado[i].X()[1] + 0.5 * step * (f[1](resultado[i]) + f[1](resultado[i + 1]));
    }

    return resultado;
};