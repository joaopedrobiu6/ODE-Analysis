#ifndef __FUNCTIONS__
#define __FUNCTIONS__

#include "pendulum_parametric.h"

class functions
{
public:
    functions() = default;
    functions(double g_,double L_);
    virtual ~functions() {};

    // Getters
    double get_g();
    double get_L();

    // Declarações de Funções
    // Draw das soluções das ODE's
    void ODE_Draw(const char *s,
                  const char *xaxis,
                  const char *yaxis,
                  int color,
                  std::vector<ODEpoint> resultado);
    // Função para o cálcul de Moving Average e função para o Draw dos gráficos das médias
    std::vector<double> Moving_Average(std::vector<ODEpoint> r_new_time, double time_step, double tw);

    void Moving_Average_Draw(const char *filename,
                             const char *xaxis,
                             const char *yaxis,
                             int color,
                             std::vector<ODEpoint> resultado,
                             std::vector<double> media,
                             double epsilon,
                             double time_step,
                             double tw);
    // Draw dos Retratos de Fase
    void Phase_Draw(const char *title,
                    const char *s,
                    int color,
                    std::vector<ODEpoint> resultado);

    // Função de verficação da estabilidade para dados delta e ep
    int Is_Stable(double *media, int n);

    // Função para determinar o mapa de estabilidade
    void Draw_Stability(const char *filename, int numero, double step);
    // Função para desenhar o mapa de estabilidade com zoom na regiao instavle+
    void Draw_Stability_Zoom(const char *filename, int numero, double step);
    // Função para desenhar o gráfico comparação entre a solução através de runge kutta e solução analitica pela análise de fourier
    void Comparacao_analitica(const char *filename, const char *xaxis, const char *yaxis, int color, std::vector<ODEpoint> resultado, std::vector<std::pair<double, double>> Freq);

    // Análise de Fourier
    std::vector<std::pair<double, double>> Fourier(std::vector<ODEpoint> resultado, double step);
    void Fourier_Draw_Harmonics(const char *s, std::vector<std::pair<double, double>>);
    void Fourier_Draw_Max(const char *s, std::vector<std::pair<double, double>>,int qual_deles);

private:
    //Constantes
    double g;                      // Aceleração Gravítica
    double L;        // Comprimento do pêndulo
};
#endif