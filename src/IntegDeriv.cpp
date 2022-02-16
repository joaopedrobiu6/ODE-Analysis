#include "IntegDeriv.h"

void IntegDeriv::firstDerivative(double x, double h, double &Derivative, std::string option)
{
    if (option == "forward")
    {
        Derivative = (-F(x + 2 * h) + 4 * F(x + h) - 3 * F(x)) / (2 * h);
    }
    if (option == "central")
    {
        Derivative = (F(x + h) - F(x - h)) / (2 * h);
    }
    if (option == "backward")
    {
        Derivative = (F(x - 2 * h) - 4 * F(x - h) + 3 * F(x)) / (2 * h);
    }
    if (option == "fivepoint")
    {
        Derivative = ((F(x - 2 * h) + 8 * F(x + h)) - (8 * F(x - h) + F(x + 2 * h))) / (12 * h);
    }
};
void IntegDeriv::secondDerivative(double x, double h, double &SndDerivative, std::string option)
{

    if (option == "threepoint")
    {
        SndDerivative = (F(x + h) - 2 * F(x) + F(x - h)) / (h * h);
    }
    if (option == "fivepoint")
    {
        SndDerivative = (-F(x - 2 * h) + 16 * F(x - h) - 30 * F(x) + 16 * F(x + h) - F(x + 2 * h)) / (12 * h * h);
    }
};
void IntegDeriv::thirdDerivative(double x, double h, double &TrdDerivative)
{

    TrdDerivative = (-F(x + 3 * h) + 8 * F(x + 2 * h) - 13 * F(x + h) + 13 * F(x - h) - 8 * F(x - 2 * h) + F(x - 3 * h)) / (8 * h * h * h);
};
void IntegDeriv::fourthDerivative(double x, double h, double &FthDerivative)
{

    FthDerivative = (-F(x - 3 * h) + 12 * F(x - 2 * h) - 39 * F(x - h) + 56 * F(x) - F(x + 3 * h) + 12 * F(x + 2 * h) - 39 * F(x + h)) / (6 * h * h * h * h);
};

/*************************************************************/
/*///////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////*/
/*///////////////////INTEGRATION METHODS////////////////////*/
/*///////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////*/
/*************************************************************/

void IntegDeriv::TrapezoidalRule(double xi, double xf, double &Integral, double &Error)
{
    ///////////////// Versão do Pedro ////////////////////////////
    int k = 0, flag = 0, i;
    double sum, n, sum2, sum3, deriv = 0, h, blank, Error2;

    while (1 && k < 20 & !flag)
    {
        k++;
        n = std::pow(2, k - 1);
        h = (xf - xi) / n;
        for (int j = 1; j <= n; j++)
        {
            secondDerivative(xi + j * h - h / 2, h / 4, sum3, "threepoint");
            deriv += sum3;
        }
        deriv /= n; // media das segundas derivadas a meio de cada fatia

        Error = std::fabs(h * h * (xf - xi) * deriv / 12); // formula do erro
        Error2 = Integral;
        if (k == 1)
        {
            Integral = (xf - xi) / 2 * (F(xi) + F(xf)); // regra do trapézio
        }
        else
        {
            sum = 0;
            sum2 = 0;
            for (i = 1; i <= n / 2; i++)
            {
                sum2 = F(xi + (2 * i - 1) * h); // formula nos slides
                if (!(sum2 > -1e5 && sum2 < 1e5))
                {
                    flag = 1;
                    break;
                }
                // std::cout << sum2 << std::endl;
                sum += sum2;
            }
            Integral = 0.5 * Integral + h * sum;
        }
        Error2 = std::fabs(Error2 - Integral); // erro2 é a diferença entre cada iteração

        if (Error < std::pow(10, -5) && Error2 < std::pow(10, -4) && Error != 0)
        { // erro nao pode ser igual a 0, seria uma coincidencia
            break;
        }
    }

    if (!flag)
    {
        std::cout << "Valor do Integral: " << Integral << " com " << k - 1 << " iterações, " << n << " partições e erros de " << Error << " e " << Error2 << std::endl;
    }
    else
    {
        std::cout << "ERRO: Existe uma singularidade perto de " << (xi + (2 * i - 1) * h) << std::endl;
    }
};

void IntegDeriv::SimpsonRule(double xi, double xf, double &Integral, double &Error)
{
    /////////////////////////////Versão do Pedro/////////////////////////////
    int k = 1, i, flag = 0;
    double deriv = 0, h, n, blank, sum, sum3, Error2;

    Integral = 0;

    while (1 && k < 20 & !flag)
    {
        k++;
        n = std::pow(2, k - 1);
        h = (xf - xi) / n;

        for (int j = 1; j <= n; j++)
        {
            fourthDerivative(xi + j * h - h / 2, h / 4, sum3);
            deriv += sum3;
        }
        deriv /= n; // média das quartas derivadas no centro de cada fatia

        Error = std::fabs(h * h * h * h * (xf - xi) * deriv / 180); // formula do erro
        Error2 = Integral;

        Integral = 0;

        for (i = 1; i <= n / 2; i++)
        {
            sum = h / 3. * (F(xi + (2 * i - 2) * h) + 4 * F(xi + (2 * i - 1) * h) + F(xi + 2 * i * h));

            if (!(sum > -1e5 && sum < 1e5))
            {
                flag = 1;
                break;
            }
            Integral += sum;
        }

        Error2 = std::fabs(Error2 - Integral);

        if (Error < std::pow(10, -5) && Error2 < std::pow(10, -4) && Error != 0)
        { // erro nao pode ser igual a 0, seria uma coincidencia
            break;
        }
    }

    if (!flag)
    {
        std::cout << "Valor do Integral: " << Integral << " com " << k - 1 << " iterações, " << n << " partições e erros de " << Error << " e " << Error2 << std::endl;
    }
    else
    {
        std::cout << "ERRO: Existe uma singularidade perto de " << (xi + (2 * i - 1) * h) << std::endl;
    }
};

void IntegDeriv::MonteCarlo_VonNeumann(double xi, double xf, double &Integral, double &Error)
{
    typedef std::chrono::high_resolution_clock Clock;
    auto t1 = Clock::now();

    double max = -1.79769e+308;
    double min = 1.79769e+308;
    double step = (xf - xi) / 1000;
    double loop = xi;
    double area;
    int num_rands = 10000000; // Numero total de pontos aleatorios gerados
    double rands_dentro_integral_positivo = 0;
    double rands_dentro_integral_negativo = 0;
    double total = 0;

    // Calcular max e min do gráf;
    while (loop <= xf)
    {
        if (F(loop) > max)
        {
            max = F(loop);
        }
        if (F(loop) < min)
        {
            min = F(loop);
        }
        loop += step;
    }

    if (min > 0)
        min = 0;
    if (max < 0)
        max = 0;

    area = (xf - xi) * std::fabs((max + step) - (min - step));

    gRandom = new TRandom3;
    gRandom->SetSeed(time(NULL));

    for (int i = 0; i < num_rands; i++)
    {
        double x = gRandom->Uniform(xi, xf);
        double y = gRandom->Uniform(min, max);
        if (y > 0 && y < F(x))
        {
            rands_dentro_integral_positivo += 1;
            total += 1;
        }
        if (y < 0 && y > F(x))
        {
            rands_dentro_integral_negativo += 1;
            total += 1;
        }
    }

    Integral = area * (double)(rands_dentro_integral_positivo / num_rands) - area * (double)(rands_dentro_integral_negativo / num_rands);
    Error = std::fabs(((xf - xi) * max) / (double)num_rands) * std::sqrt(total * (1 - (double)(total / num_rands)));

    auto t2 = Clock::now();
    std::cout << "Tempo de Execucao: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " milisegundos" << std::endl;

    std::cout << "Numero total de pontos gerados: " << num_rands << std::endl;
    std::cout << "Area total positiva: " << area * (double)(rands_dentro_integral_positivo / num_rands) << std::endl;
    std::cout << "Area total negativa: " << area * (double)(rands_dentro_integral_negativo / num_rands) << std::endl;
};

void IntegDeriv::MonteCarlo_VonNeumann_Trapezoidal(double xi, double xf, double &Integral, double &Error)
{
    typedef std::chrono::high_resolution_clock Clock;
    auto t1 = Clock::now();

    double step = (xf - xi) / 1000;
    int num_rands = 10000000; // Numero total de pontos aleatorios gerados
    int num_particoes = 100;
    double step_particoes = (xf - xi) / num_particoes;
    double total = 0;
    double area_total = 0;

    Integral = 0;
    Error = 0;

    for (int k = 0; k < num_particoes; k++)
    {
        double max = -1.79769e+308;
        double min = 1.79769e+308;
        double loop = xi + k * step_particoes;
        double area;
        double rands_dentro_integral_positivo = 0;
        double rands_dentro_integral_negativo = 0;

        // Calcular max e min do gráf, dentro de cada particao;
        while (loop <= xi + (k + 1) * step_particoes)
        {
            if (F(loop) > max)
            {
                max = F(loop);
            }
            if (F(loop) < min)
            {
                min = F(loop);
            }
            loop += step;
        }

        if (min > 0)
            min = 0;
        if (max < 0)
            max = 0;

        area = ((xi + (k + 1) * step_particoes) - (xi + k * step_particoes)) * std::fabs((max + step) - (min - step));

        area_total += area;

        gRandom = new TRandom3;
        gRandom->SetSeed(time(NULL));

        for (int i = 0; i < num_rands / num_particoes; i++)
        {
            double x = gRandom->Uniform(xi + k * step_particoes, xi + (k + 1) * step_particoes);
            double y = gRandom->Uniform(min, max);

            if (y > 0 && y < F(x))
            {
                rands_dentro_integral_positivo += 1;
                total += 1;
            }
            if (y < 0 && y > F(x))
            {
                rands_dentro_integral_negativo += 1;
                total += 1;
            }
        }

        Integral = Integral + area * (double)(rands_dentro_integral_positivo / (num_rands / num_particoes)) - area * (double)(rands_dentro_integral_negativo / (num_rands / num_particoes));
    }

    Error = Error + std::fabs(area_total / (double)num_rands) * std::sqrt(total * (1 - (double)(total / num_rands)));

    auto t2 = Clock::now();
    std::cout << "Tempo de Execucao: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " milisegundos" << std::endl;

    std::cout << "Numero total de pontos gerados: " << num_rands << std::endl;
    std::cout << "Numero total de particoes: " << num_particoes << std::endl;
};
