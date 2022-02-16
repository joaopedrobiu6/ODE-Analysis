#ifndef __INTEGDERIV__
#define __INTEGDERIV__
#include <iomanip>
#include <random>
#include <chrono>
#include "TRandom3.h"

#include "Functor.h"

class IntegDeriv
{
public:
    IntegDeriv(Functor &Funk) : F(Funk) { ; }
    ~IntegDeriv() = default;
    // integration methods
    void TrapezoidalRule(double xi, double xf, double &Integral, double &Error);
    void SimpsonRule(double xi, double xf, double &Integral, double &Error);
    void MonteCarlo_VonNeumann(double xi, double xf, double &Integral, double &Error);
    void MonteCarlo_VonNeumann_Trapezoidal(double xi, double xf, double &Integral, double &Error);

    // derivative methods
    void firstDerivative(double x, double h, double &Derivative, std::string option);
    void secondDerivative(double x, double h, double &SndDerivative, std::string option);
    void thirdDerivative(double x, double h, double &TrdDerivative);
    void fourthDerivative(double x, double h, double &FthDerivative);

private:
    Functor &F;
};
#endif