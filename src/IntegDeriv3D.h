#ifndef __INTEGDERIV3D__
#define __INTEGDERIV3D__

#include "TRandom3.h"

#include <iostream>
#include <math.h>

class IntegDeriv3D
{
public: 
    //Construtor
    IntegDeriv3D() = default;
    //Destructor
    ~IntegDeriv3D() = default;

    //Overload do operator() para definir a formula para o cálculo da Irradiância
    double operator()(double x, double y);

    double fourthDerivative1(double x, double y, double h);
    double fourthDerivative2(double xi, double xf, double y, double h);

    //Funções para calcular os integrais, retornam o valor do integral e o erro
    double SimpsonRule1(double xi, double xf, double y);
    std::pair<double, double> SimpsonRule2(double xi, double xf, double yi, double yf);
    std::pair<double, double> MonteCarlo();
};

#endif
