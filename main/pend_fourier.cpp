#include "functions.h"

int main()
{
    std::cout << "PROJETO GRUPO A03\nGuilherme Lourinho 100310\nJoão Biu 100335\nMiguel Lameiras 100354\nPedro Fialho 100364\n"
              << std::endl;

    functions funcoes(9.8, 1);

    // Cálculo da Solução para (delta, epsilon) = (0.7, 0.2)
    double delta = 0.7,
           epsilon = 0.2;
    double time_step = 1e-3;
    double tw_real = 8;
    double tw = tw_real / time_step;

    std::cout << "\n***** Alínea b1) *****" << std::endl;
    std::cout << "Resolvendo a Equação Diferencial com o Método Runge Kutta de Quarta Ordem..." << std::endl;
    // Condições iniciais: posição inicial 0.1 - velocidade inicial 0
    pendulum_parametric pendulum(funcoes.get_L(), {0.1, 0});

    // Input das Equações do Movimento do Pêndulo (ODE, já adimensionais)
    pendulum.SetFunction(0, [](ODEpoint p)
                         { return p.X()[1]; });

    pendulum.SetFunction(1, [&](ODEpoint p)
                         { return -((delta + epsilon * std::cos(p.T())) * std::sin((p.X())[0])); });

    // Cálculo da Solução
    std::vector<ODEpoint> resultado = pendulum.RungeKutta4(100, time_step);

    std::cout << "\n***** Alínea b2) *****" << std::endl;
    // Draw do gráfico das soluções
    funcoes.ODE_Draw("graph1.png", "time", "#hat{#theta}", 50, resultado);

    // Dimensionalizar o tempo para o grafico de Theta(t), [t] = segundos
    std::vector<ODEpoint> r_new_time = resultado;
    for (int i = 0; i < resultado.size(); i++)
    {
        r_new_time[i].T() *= 1 / sqrt(funcoes.get_g() / (funcoes.get_L() * delta));
    }

    // Draw do novo gráfico com o tempo em segundos
    funcoes.ODE_Draw("graph2.png", "time [s]", "#theta [rad]", 62, r_new_time);

    std::cout << "\n***** Alínea b3) *****" << std::endl;

    double amp_max = 0, amp_min = DBL_MAX; // Variáveis que guardam as amplitudes maximas e minimas do grafico
    std::vector<double> amp_all;

    for (int i = 1; i < r_new_time.size() - 1; i++)
    {
        if (fabs(r_new_time[i - 1].X()[0]) <= fabs(r_new_time[i].X()[0]) && fabs(r_new_time[i].X()[0]) >= fabs(r_new_time[i + 1].X()[0]))
        {
            amp_all.push_back(fabs(r_new_time[i].X()[0]));

            if (fabs(r_new_time[i].X()[0]) > amp_max)
            {
                amp_max = fabs(r_new_time[i].X()[0]);
            }
            else if (fabs(r_new_time[i].X()[0]) < amp_min)
            {
                amp_min = fabs(r_new_time[i].X()[0]);
            }
        }
    }

    // Cálculo da Profundidade de Modulação e prints
    double profundidade_modulacao = (amp_max - amp_min) / (amp_max + amp_min);

    std::cout << "\nAmplitude máxima: " << amp_max << std::endl;
    std::cout << "\nAmplitude mínima: " << amp_min << std::endl;
    std::cout << "\nProfundidade de Modulação: " << profundidade_modulacao << std::endl;

    std::cout << "\n***** Alínea c1) *****" << std::endl;
    // MOVING AVERAGES
    // Primeiros valores de delta e epsilon
    delta = 0.4;
    epsilon = 0.2;

    pendulum_parametric pendulum1(funcoes.get_L(), {0.1, 0}); // criar novo objeto para resolver com os novos parâmetros

    // Input das equações com os novos parâmetros
    pendulum1.SetFunction(0, [](ODEpoint p)
                          { return p.X()[1]; });

    pendulum1.SetFunction(1, [&](ODEpoint p)
                          { return -((delta + epsilon * std::cos(p.T())) * std::sin((p.X())[0])); });

    // Resolver a Equação Diferencial com RK4 na classe "pendulum_parametric"
    std::vector<ODEpoint> resultado1 = pendulum1.RungeKutta4(150, time_step);

    // Cálculo dos valores médios com o método Moving Average
    std::vector<double> mov_avg = funcoes.Moving_Average(resultado1, time_step, tw);
    // Draw do grafico dos valores medios de hat{theta}
    funcoes.Moving_Average_Draw("avg_draw1.png", "time", "#hat{#theta}", 152, resultado1, mov_avg, epsilon, time_step, tw);
    // Diagrama de fase para (delta, epsilon) = (0.4, 0.2)
    funcoes.Phase_Draw("Retrato de Fase para #delta = 0.4 #varepsilon = 0.2", "phase1.png", 50, resultado1);

    // Segundos valores de delta e epsilon
    delta = 0.4;
    epsilon = 0.36;

    pendulum_parametric pendulum2(funcoes.get_L(), {0.1, 0}); // Novo objeto para resolver as eqs diferenciais com os parametros (delta, epsilon) = (0.4, 0.36)

    // Input das equações com os novos parâmetros
    pendulum2.SetFunction(0, [](ODEpoint p)
                          { return p.X()[1]; });

    pendulum2.SetFunction(1, [&](ODEpoint p)
                          { return -((delta + epsilon * std::cos(p.T())) * std::sin((p.X())[0])); });

    // Resolver a Equação Diferencial com RK4 na classe "pendulum_parametric"
    std::vector<ODEpoint> resultado2 = pendulum2.RungeKutta4(150, time_step);

    // Cálculo dos valores médios com o método Moving Average
    std::vector<double> mov_avg2 = funcoes.Moving_Average(resultado2, time_step, tw);
    // Draw do grafico dos valores medios de hat{theta}
    funcoes.Moving_Average_Draw("avg_draw2.png", "time", "#hat{#theta}", 152, resultado2, mov_avg2, epsilon, time_step, tw);
    // Diagrama de fase para (delta, epsilon) = (0.4, 0.36)
    funcoes.Phase_Draw("Retrato de Fase para #delta = 0.4 #varepsilon = 0.36", "phase2.png", 62, resultado2);

    std::cout << "\n***** Alínea c2) *****" << std::endl;
    // No relatorio foram gerados muito mais pontos para aumentar a resolução, no entanto, devido ao extenso tempo de computação, o número foi reduzido
    // Diagrama de pontos da estabilidade da solução da eq diferencial para os diferentes valores de delta e epsilon
    funcoes.Draw_Stability("stability_draw4.png", 350, 0.001);
    // Função para imprimir um diagrama com limites mais pequenos (dar zoom ao grafico anterior para complementar o relatório)
    //funcoes.Draw_Stability_Zoom("stability_map_zoom.png", 350, 0.05);

    std::cout << "\n***** Alínea d) *****" << std::endl;
    auto Xk = funcoes.Fourier(r_new_time, time_step * sqrt(0.7 / funcoes.get_g()));
    funcoes.Fourier_Draw_Harmonics("fourier1.png", Xk);
    funcoes.Fourier_Draw_Max("fourier1_Max.png", Xk, 90);
    funcoes.Comparacao_analitica("analitica.png", "time", "#hat{#theta}", 50, r_new_time, Xk);

    std::vector<ODEpoint> r_new_time2 = resultado2;
    for (int i = 0; i < resultado2.size(); i++)
    {
        r_new_time2[i].T() *= 1 / sqrt(funcoes.get_g() / (funcoes.get_L() * delta));
    }

    Xk = funcoes.Fourier(r_new_time2, time_step * sqrt(0.4 / funcoes.get_g()));
    funcoes.Fourier_Draw_Max("fourier2_Max.png", Xk, 0);

    return 0;
}
