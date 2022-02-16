#include "IntegDeriv3D.h"
/*
    Este operador() é a função da irradiância que vai ser integrada em todo o plano, para obter a potência incidente.
    É uma função de duas variáveis, x e y, visto que z está fixo.
*/
double IntegDeriv3D::operator()(double x, double y)
{
    double fluxo = 5e3;
    double z = 100;

    return sin(x)*sin(x)*cos(y);
};

/*
    Este é o cálculo do integral da função da irradiância num dado y fixo, em todo o domínio de x. É usado o método de Simpson.
    É ainda feita uma abordagem iterativa, calculando o integral várias vezes em divisões do domínio de x em potências de 2, até o Erro ser inferior a 1e-5.
    O Erro é obtido através da quarta derivada da função a integrar. É ainda usado um Erro secundário Error2 para ter a certeza que os integrais estão a convergir.
    Os erros não são retornados mas influenciam o cálculo do Integral.
*/

double IntegDeriv3D::SimpsonRule1(double xi, double xf, double y)
{

    int k = 1, i, flag = 0;
    double deriv = 0, h, n, blank, sum, Error2;
    double Integral = 0, Error;

    while (1 && k < 20 & !flag)                          //início do ciclo das divisões do domínio
    {
        k++;
        n = std::pow(2, k - 1);                         //as divisões do domínio são potências de 2
        h = (xf - xi) / n;

        for (int j = 1; j <= n; j++)
        {
            deriv += fourthDerivative1(xi + j * h - h / 2, y, h / 4);
        }
        deriv /= n;                                    // média das quartas derivadas no centro de cada "fatia"

        Error = std::fabs(h * h * h * h * (xf - xi) * deriv / 180);     //formula do erro, com a quarta derivada
        Error2 = Integral;

        Integral = 0;

        for (i = 1; i <= n / 2; i++)
        {
            sum = h / 3. * (operator()(xi + (2 * i - 2) * h, y) + 4 * operator()(xi + (2 * i - 1) * h, y) + operator()(xi + 2 * i * h, y));   //formula da regra de Simpson

            if (!(sum > -1e6 && sum < 1e6))
            {
                flag = 1;                                  //uma flag alerta caso a função esteja perto de uma singularidade
                break;
            }
            Integral += sum;
        }

        Error2 = std::fabs(Error2 - Integral);

        if (Error < std::pow(10, -5) && Error2 < std::pow(10, -4) && Error != 0)  //Fecho do ciclo quando o Erro é inferior a 10-5
        { 
            break;
        }
    }

    if (flag)
    {
        std::cout << "ERRO: Existe uma singularidade perto de (" << (xi + (2 * i - 1) * h) << ", " << y << ")" << std::endl;
        return -1;
    }


    return Integral;
};


/*
    Este é o cálculo do integral em y dos integrais já obtidos em x. É usado o método de Simpson como anteriormente, mas a função a integrar é ela mesmo um integral de Simpson.
    O Erro é obtido através de uma quarta derivada dessa função. A função retorna o valor do integral e o erro.
*/

std::pair<double, double> IntegDeriv3D::SimpsonRule2(double xi, double xf, double yi, double yf)
{
    /////////////////////////////Versão do Pedro/////////////////////////////
    int k = 1, i, flag = 0;
    double deriv = 0, h, n, blank, sum, sum3, Error2;

    double Integral, Error;

    Integral = 0;

    while (1 && k < 20 & !flag)                               //início do ciclo das divisões do domínio
    {
        k++;
        n = std::pow(2, k - 1);                             //as divisões do domínio são potências de 2
        h = (xf - xi) / n;

        Error = 0;

        for (int j = 1; j <= n; j++)
        {
            sum3 = fourthDerivative2(xi, xf, yi + j * h - h / 2, h / 4);
            deriv += sum3;
        }
        deriv /= n;                                        // média das quartas derivadas no centro de cada "fatia"

        Error = std::fabs(h * h * h * h * (yf - yi) * deriv / 180);     //formula do erro, com a quarta derivada
        Error2 = Integral;

        Integral = 0;

        for (i = 1; i <= n / 2; i++)
        {
            sum = h / 3. * (SimpsonRule1(xi, xf, yi + (2 * i - 2) * h) + 4 * SimpsonRule1(xi, xf, yi + (2 * i - 1) * h) + SimpsonRule1(xi, xf, yi + 2 * i * h));

            if (!(sum > -1e6 && sum < 1e6))
            {
                flag = 1;
                break;
            }
            Integral += sum;
        }

        Error2 = std::fabs(Error2 - Integral);

        if (Error < std::pow(10, -3) && Error2 < std::pow(10, -2) && Error != 0)            //Fecho do ciclo quando o Erro é inferior a 10-5
        { 
            break;
        }
    }

    if (flag)
    {
        std::cout << "ERRO: Existe uma singularidade perto de y= " << yi + (2 * i - 1) * h << std::endl;
        return std::make_pair(-1., -1.);
    }

    return std::make_pair(Integral, Error);
};



/*
    Este é o cálculo da quarta derivada parcial na direção de x da função da irradiância
*/
double IntegDeriv3D::fourthDerivative1(double x, double y, double h)
{
    return (-operator()(x - 3 * h, y) + 12 * operator()(x - 2 * h, y) - 39 * operator()(x - h, y) + 56 * operator()(x, y) - operator()(x + 3 * h, y) + 12 * operator()(x + 2 * h, y) - 39 * operator()(x + h, y)) / (6 * h * h * h * h);
};


/*
    Este é o cálculo da quarta derivada em y da função SimpsonRule1, que é por sua vez o cálculo de um integral
*/
double IntegDeriv3D::fourthDerivative2(double xi, double xf, double y, double h)
{
    return (-SimpsonRule1(xi, xf, y - 3 * h) + 12 * SimpsonRule1(xi, xf, y - 2 * h) - 39 * SimpsonRule1(xi, xf, y - h) + 56 * SimpsonRule1(xi, xf, y) - SimpsonRule1(xi, xf, y + 3 * h) + 12 * SimpsonRule1(xi, xf, y + 2 * h) - 39 * SimpsonRule1(xi, xf, y + h)) / (6 * h * h * h * h);
};

std::pair<double, double> IntegDeriv3D::MonteCarlo()
{
    double num_rands = 500000000; // Número de Pontos Aleatórios gerados
    double max = 0; //Máximo da função para que se possa definir o tamanho da região
    double num_rands_dentro = 0; //Número de pontos aleatórios que estão abaixo da curva da função
    std::pair<double, double> resultado; //par de valores, o primeiro o integral, o segundo o erro

    //Calcular o máximo da função
    for (int i = -200; i < 200; i++)
    {
        for (int j = -200; j < 200; j++)
        {
            if (operator()(i, j) > max)
            {
                max = operator()(i, j);
            }
        }
    }

    //Volume total da região onde vão ser gerados os pontos aleatórios
    double volume = 400 * 400 * max;

    gRandom = new TRandom3;
    gRandom->SetSeed(time(NULL));

    // Método monte carlo
    for (int i = 0; i < num_rands; i++)
    {
        //Gerar pontos aleatorios
        double randx = gRandom->Uniform(-200,200);
        double randy = gRandom->Uniform(-200,200);
        double randz = gRandom->Uniform(0,max);
        //Verificar se o ponto gerado está abaixo da curva da função
        if (randz <= operator()(randx, randy))
        {
            num_rands_dentro++;
        }
    }

    resultado.first = volume * (num_rands_dentro/num_rands); //Valor do Integral pelo método Monte Carlo Von Neumman
    resultado.second = (resultado.first/num_rands_dentro); //Expressão para o cálculo do erro
    return resultado;
};