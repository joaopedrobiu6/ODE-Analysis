#include "pendulum.h"

pendulum::pendulum(double l, const Xvar &X)
{
    X0 = X;
    L = l;
    f[0] = [](ODEpoint p) { // definir a lambda function
        auto temp = p.X();
        return temp[1]; /*return velocidade*/
    };
    f[1] = [&](ODEpoint p)
    {
        auto temp = p.X();
        return (((-9.8 / L) + std::cos(p.T())) * std::sin((p.X()[0]))); /*return aceleração*/
    };
};

pendulum::pendulum(double l, const std::initializer_list<double> &X)
{
    Xvar var(X);
    X0 = var;
    L = l;
    f[0] = [](ODEpoint p)
    {
        auto temp = p.X();
        return temp[1]; /*return velocidade*/
    };
    f[1] = [&](ODEpoint p)
    {
        auto temp = p.X();
        return (((-9.8 / L) + std::cos(p.T())) * std::sin((p.X())[0])); /*return aceleração*/
    };
};

pendulum::pendulum(double l, double theta_0, double theta_vel_0)
{
    std::vector<double> temp;
    temp[0] = theta_0;
    temp[1] = theta_vel_0;
    temp[2] = 0;
    Xvar var(temp);
    X0 = var;
    L = l;
    f[0] = [](ODEpoint p)
    {
        auto temp = p.X();
        return temp[1]; /*return velocidade*/
    };
    f[1] = [&](ODEpoint p)
    {
        auto temp = p.X();
        return (((-9.8 / L) + std::cos(p.T())) * std::sin((p.X())[0])); /*return aceleração*/
    };
};

const std::vector<ODEpoint> &pendulum::VerletSolver(double Time, double step)
{
    std::vector<ODEpoint> info3;
    int n = Time / step;

    // resizes
    info3.resize(n);
    int n1 = X0.X().size(); // n1 é o tamanho de x de X0 - número de variáveis dependentes

    for (int i = 0; i < n; i++)
        info3[i].X().resize(n1);

    // condições iniciais
    info3[0].T() = 0;

    for (int i = 0; i < n1; i++)
        info3[0].X()[i] = X0.X()[i];

    // primeira iteração
    info3[1].T() = step + info3[0].T();
    info3[1].X()[0] = info3[0].X()[0] + step * f[0](info3[0]);
    info3[1].X()[1] = info3[0].X()[1] + step * f[1](info3[0]);

    for (int i = 1; i < n - 1; i++)
    {
        info3[i + 1].T() = step + info3[i].T();
        info3[i + 1].X()[0] = info3[i - 1].X()[0] + 2 * step * f[0](info3[i]);
        info3[i + 1].X()[1] = info3[i - 1].X()[1] + 2 * step * f[1](info3[i]);
    }

    MS.insert({"verlet", info3});
    return MS["verlet"];
};

const std::vector<ODEpoint> &pendulum::LeapFrogImprovedSolver(double Time, double step)
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
    info4[1].X()[0] = info4[0].X()[0] + step * f[0](info4[0]);
    info4[1].X()[1] = info4[0].X()[1] + step * f[1](info4[0]);

    for (int i = 1; i < n - 1; i++)
    {
        info4[i + 1].T() = step + info4[i].T();
        info4[i + 1].X()[0] = info4[i].X()[0] + step * f[0](info4[i]) + 0.5 * step * step * f[1](info4[i]);
        info4[i + 1].X()[1] = info4[i].X()[1] + 0.5 * step * (f[1](info4[i]) + f[1](info4[i + 1]));
    }

    MS.insert({"leapfrog", info4});
    return MS["leapfrog"];
};

const std::vector<ODEpoint> &pendulum::EulerSolver(double Time, double step)
{
    std::vector<ODEpoint> info1;
    int n = Time / step;

    // resizes
    info1.resize(n);
    int n1 = X0.X().size(); // n1 é o tamanho de x de X0 - número de variáveis dependentes

    for (int i = 0; i < n; i++)
        info1[i].X().resize(n1);

    // condições iniciais
    info1[0].T() = 0;

    for (int i = 0; i < n1; i++)
        info1[0].X()[i] = X0.X()[i];

    // calculo das soluções

    for (int i = 0; i < n - 1; i++)
    {
        info1[i + 1].T() = step + info1[i].T();
        info1[i + 1].X()[0] = info1[i].X()[0] + step * f[0](info1[i]);
        info1[i + 1].X()[1] = info1[i].X()[1] + step * f[1](info1[i]);
    }

    MS.insert({"euler", info1});
    return MS["euler"];
}

const std::vector<ODEpoint> &pendulum::TrapezoidalSolver(double Time, double step)
{
    std::vector<ODEpoint> info2;
    int n = Time / step;

    // resizes
    info2.resize(n);
    int n1 = X0.X().size(); // n1 é o tamanho de x de X0 - número de variáveis dependentes

    for (int i = 0; i < n; i++)
        info2[i].X().resize(n1);

    // condições iniciais
    info2[0].T() = 0;

    for (int i = 0; i < n1; i++)
        info2[0].X()[i] = X0.X()[i];

    // calculo das soluções

    for (int i = 0; i < n - 1; i++)
    {
        info2[i + 1].T() = step + info2[i].T();
        info2[i + 1].X()[1] = info2[i].X()[1] + step * f[1](info2[i]);
        info2[i + 1].X()[0] = info2[i].X()[0] + step * 0.5 * (info2[i + 1].X()[1] + info2[i].X()[1]);
    }

    MS.insert({"trapezoidal", info2});
    return MS["trapezoidal"];
};

const std::vector<ODEpoint> &pendulum::EulerCromerSolver(double Time, double step)
{
    std::vector<ODEpoint> info5;
    int n = Time / step;

    // resizes
    info5.resize(n);
    int n1 = X0.X().size(); // n1 é o tamanho de x de X0 - número de variáveis dependentes

    for (int i = 0; i < n; i++)
        info5[i].X().resize(n1);

    // condições iniciais
    info5[0].T() = 0;

    for (int i = 0; i < n1; i++)
        info5[0].X()[i] = X0.X()[i];

    // calculo das soluções

    for (int i = 0; i < n - 1; i++)
    {
        info5[i + 1].T() = step + info5[i].T();
        info5[i + 1].X()[1] = info5[i].X()[1] + step * f[1](info5[i]);      // velocidade
        info5[i + 1].X()[0] = info5[i].X()[0] + step * info5[i + 1].X()[1]; // posição
    }

    MS.insert({"eulercromer", info5});
    return MS["eulercromer"];
};

const std::vector<ODEpoint> &pendulum::RungeKutta2(double Time, double step)
{
    std::vector<ODEpoint> info6;
    std::vector<double> K1, K2;

    // resizes
    int n = Time / step;
    int n1 = X0.X().size();

    info6.resize(n);
    K1.resize(n1);
    K2.resize(n1);

    for (int i = 0; i < n; i++)
    {
        info6[i].X().resize(n1);
    }

    // condições iniciais
    info6[0].T() = 0;

    for (int i = 0; i < n1; i++)
    {
        info6[0].X()[i] = X0.X()[i];
    }

    ODEpoint p0;
    Xvar xvar_temp;
    double tempo;
    p0.X().resize(n1);
    xvar_temp.X().resize(n1);

    for (int i = 0; i < n - 1; i++)
    {
        info6[i + 1].T() = step + info6[i].T();

        // Cálculo dos K1
        K1[0] = step * f[0](info6[i]);
        K1[1] = step * f[1](info6[i]);

        /* Cálculo dos K2
        Como K2 = step*f(ti + step/2, yi + K1/2)
        temos de criar um ODEpoint (t, y) com as coordenadas (ti + step/2, yi + K1/2)            muito verdade!!
        (dois na verdade, para as duas funções)
        */

        xvar_temp.X()[0] = info6[i].X()[0] + (0.5 * K1[0]);
        xvar_temp.X()[1] = info6[i].X()[1] + (0.5 * K1[1]);
        tempo = info6[i].T() + (step * 0.5);

        p0.SetODEpoint(tempo, xvar_temp);

        K2[0] = step * f[0](p0);
        K2[1] = step * f[1](p0);

        // Cálculo das soluções
        info6[i + 1].X()[0] = info6[i].X()[0] + K2[0];
        info6[i + 1].X()[1] = info6[i].X()[1] + K2[1];

        // std::cout << "ola3 " << i << std::endl;
        K2.clear();
        K1.clear();
    }

    MS.insert({"rungekutta2", info6});
    return MS["rungekutta2"];
};

const std::vector<ODEpoint> &pendulum::RungeKutta4(double Time, double step)
{
    std::vector<ODEpoint> info7;
    std::vector<double> K1, K2, K3, K4;

    // resizes
    int n = Time / step;
    int n1 = X0.X().size();

    info7.resize(n);
    K1.resize(n1);
    K2.resize(n1);
    K3.resize(n1);
    K4.resize(n1);

    for (int i = 0; i < n; i++)
    {
        info7[i].X().resize(n1);
    }

    // condições iniciais
    info7[0].T() = 0;

    for (int i = 0; i < n1; i++)
    {
        info7[0].X()[i] = X0.X()[i];
    }

    ODEpoint p0;
    Xvar xvar_temp;
    double tempo;
    p0.X().resize(n1);
    xvar_temp.X().resize(n1);

    for (int i = 0; i < n - 1; i++)
    {
        info7[i + 1].T() = info7[i].T() + step;

        // Cálculo dos K1
        K1[0] = step * f[0](info7[i]);
        K1[1] = step * f[1](info7[i]);

        // Cálculo dos K2
        xvar_temp.X()[0] = info7[i].X()[0] + (0.5 * step * K1[0]);
        xvar_temp.X()[1] = info7[i].X()[1] + (0.5 * step * K1[1]);
        tempo = info7[i].T() + (step * 0.5);

        p0.SetODEpoint(tempo, xvar_temp);

        K2[0] = 2 * step * f[0](p0);
        K2[1] = 2 * step * f[1](p0);

        // Cálculo dos K3
        xvar_temp.X()[0] = info7[i].X()[0] + (0.5 * step * K2[0]);
        xvar_temp.X()[1] = info7[i].X()[1] + (0.5 * step * K2[1]);
        tempo = info7[i].T() + (step * 0.5);

        p0.SetODEpoint(tempo, xvar_temp);

        K3[0] = 2 * step * f[0](p0);
        K3[1] = 2 * step * f[1](p0);

        // Cálculo dos K4
        xvar_temp.X()[0] = info7[i].X()[0] + (step * K2[0]);
        xvar_temp.X()[1] = info7[i].X()[1] + (step * K2[1]);
        tempo = info7[i].T() + (step);

        p0.SetODEpoint(tempo, xvar_temp);

        K4[0] = step * f[0](p0);
        K4[1] = step * f[1](p0);

        // Cálculo das soluções
        info7[i + 1].X()[0] = info7[i].X()[0] + (K1[0] + K2[0] + K3[0] + K4[0]) / 6;
        info7[i + 1].X()[1] = info7[i].X()[1] + (K1[1] + K2[1] + K3[1] + K4[1]) / 6;

        K1.clear();
        K2.clear();
        K3.clear();
        K4.clear();
    }

    MS.insert({"rungekutta4", info7});
    return MS["rungekutta4"];
};

void pendulum::Draw(std::string s, double ti, double tf)
{

    TApplication *app = new TApplication("app", nullptr, nullptr);
    TCanvas *c = new TCanvas("canvas", "Pendulum", 0, 0, 1280, 720);

    TRootCanvas *r = (TRootCanvas *)c->GetCanvasImp();
    r->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

    int flag = 0;

    if (s == "euler")
    {
        this->EulerSolver(tf - ti, 1e-3);
    }
    else if (s == "trapezoidal")
    {
        this->TrapezoidalSolver(tf - ti, 1e-3);
    }
    else if (s == "verlet")
    {
        this->VerletSolver(tf - ti, 1e-3);
    }
    else if (s == "leapfrog")
    {
        this->LeapFrogImprovedSolver(tf - ti, 1e-3);
    }
    else if (s == "eulercromer")
    {
        this->EulerCromerSolver(tf - ti, 1e-3);
    }
    else if (s == "rungekutta2")
    {
        this->RungeKutta2(tf - ti, 1e-3);
    }
    else if (s == "rungekutta4")
    {
        this->RungeKutta4(tf - ti, 1e-3);
    }
    else
    {
        flag = 1;

        TMultiGraph *mg = new TMultiGraph();
        TGraph *gr[7];

        this->EulerSolver(tf - ti, 1e-3);
        this->TrapezoidalSolver(tf - ti, 1e-3);
        this->VerletSolver(tf - ti, 1e-3);
        this->LeapFrogImprovedSolver(tf - ti, 1e-3);
        this->EulerCromerSolver(tf - ti, 1e-3);
        this->RungeKutta2(tf - ti, 1e-3);
        this->RungeKutta4(tf - ti, 1e-3);

        int n = MS["euler"].size();

        double tempo[n];
        double posicao[n];
        double velocidade[n];

        for (int i = 0; i < n; i++)
        {
            tempo[i] = MS["euler"][i].T();
            posicao[i] = MS["euler"][i].X()[0];
            velocidade[i] = MS["euler"][i].X()[1];
        }
        gr[0] = new TGraph(n, tempo, posicao);

        for (int i = 0; i < n; i++)
        {
            tempo[i] = MS["trapezoidal"][i].T();
            posicao[i] = MS["trapezoidal"][i].X()[0];
            velocidade[i] = MS["trapezoidal"][i].X()[1];
        }
        gr[1] = new TGraph(n, tempo, posicao);

        for (int i = 0; i < n; i++)
        {
            tempo[i] = MS["verlet"][i].T();
            posicao[i] = MS["verlet"][i].X()[0];
            velocidade[i] = MS["verlet"][i].X()[1];
        }
        gr[2] = new TGraph(n, tempo, posicao);

        for (int i = 0; i < n; i++)
        {
            tempo[i] = MS["leapfrog"][i].T();
            posicao[i] = MS["leapfrog"][i].X()[0];
            velocidade[i] = MS["leapfrog"][i].X()[1];
        }
        gr[3] = new TGraph(n, tempo, posicao);

        for (int i = 0; i < n; i++)
        {
            tempo[i] = MS["eulercromer"][i].T();
            posicao[i] = MS["eulercromer"][i].X()[0];
            velocidade[i] = MS["eulercromer"][i].X()[1];
        }
        gr[4] = new TGraph(n, tempo, posicao);

        for (int i = 0; i < n; i++)
        {
            tempo[i] = MS["rungekutta2"][i].T();
            posicao[i] = MS["rungekutta2"][i].X()[0];
            velocidade[i] = MS["rungekutta2"][i].X()[1];
        }
        gr[5] = new TGraph(n, tempo, posicao);

        for (int i = 0; i < n; i++)
        {
            tempo[i] = MS["rungekutta4"][i].T();
            posicao[i] = MS["rungekutta4"][i].X()[0];
            velocidade[i] = MS["rungekutta4"][i].X()[1];
        }
        gr[6] = new TGraph(n, tempo, posicao);

        // mg->Add(gr[0], "cp");
        // mg->Add(gr[1], "cp");
        // mg->Add(gr[2], "cp");
        // mg->Add(gr[3], "cp");
        // mg->Add(gr[4], "cp");
        mg->Add(gr[5], "cp");
        mg->Add(gr[6], "cp");

        gROOT->SetStyle("Plain");

        gr[0]->SetLineColor(80);
        gr[0]->SetLineWidth(2);
        gr[0]->SetTitle("Euler");

        gr[1]->SetLineColor(20);
        gr[1]->SetLineWidth(2);
        gr[1]->SetTitle("Trapezoidal");

        gr[2]->SetLineColor(30);
        gr[2]->SetLineWidth(2);
        gr[2]->SetTitle("Verlet");

        gr[3]->SetLineColor(40);
        gr[3]->SetLineWidth(2);
        gr[3]->SetTitle("Leapfrog");

        gr[4]->SetLineColor(50);
        gr[4]->SetLineWidth(2);
        gr[4]->SetTitle("eulercromer");

        gr[5]->SetLineColor(60);
        gr[5]->SetLineWidth(2);
        gr[5]->SetTitle("rungekutta2");

        gr[6]->SetLineColor(70);
        gr[6]->SetLineWidth(2);
        gr[6]->SetTitle("rungekutta4");

        mg->SetTitle("Pendulo;Tempo;Posicao");
        gr[6]->Draw("ALP");
        mg->Draw("LP");

        c->Update();
        c->BuildLegend();
        app->Run();
    }

    if (!flag)
    {
        int n = MS[s].size();

        double tempo[n];
        double posicao[n];
        double velocidade[n];

        for (int i = 0; i < n; i++)
        {
            tempo[i] = MS[s][i].T();
            posicao[i] = MS[s][i].X()[0];
            // std::cout << "\nposição para i: " << posicao[i] << std::endl;
            velocidade[i] = MS[s][i].X()[1];
        }

        TGraph *gr = new TGraph(n, tempo, posicao);
        gROOT->SetStyle("Plain");
        gr->SetLineColor(61);
        gr->SetLineWidth(4);

        gr->SetTitle("Pendulo;Tempo;Posicao");
        gr->Draw("AC");

        c->Update();
        app->Run();
    }

    std::cout << "\nlua lula" << std::endl;
};