#include "MonteCarlo.h"

void MonteCarlo::MCIntegration(double xi, double xf, double &Integral, double &Error)
{
    double sum = 0., negsum = 0., err = 0., var = 0.;
    double max = -DBL_MAX, min = DBL_MAX;
    double step = 1E-3, search = xi;
    gRandom = new TRandom3(0);

    // Get function Max
    while (search <= xf)
    {
        if (F(search) > max)
            max = F(search);
        if (F(search) < min)
            min = F(search);
        search = search + step;
    }

    // Get Area
    double area = (xf - xi) * std::fabs(max - min);

    for (int i = 0; i < N; i++)
    {
        double x_ = gRandom->Uniform(xi, xf);
        double y_ = gRandom->Uniform(min, max);
        if (y_ > 0 && y_ < F(x_))
            sum = sum + 1;
        if (y_ < 0 && y_ > F(x_))
            negsum = negsum + 1;
    }
    Integral = area * (double)(sum / N) - area * (double)(negsum / N);
    Error = (area / N) * sqrt((sum - negsum) * (1 - ((sum - negsum) / N)));
};