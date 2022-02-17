#include "Fitter.h"

int main()
{

  std::vector<std::vector<float>> vec = {{0.425, 1.826}, {1.489, 2.993}, {2.979, 1.324}};

  TF1 *f_1 = new TF1("f", "[0]*sin(x)+[1]", 0, 10);
  f_1->SetParNames("a", "b");
  f_1->SetParameters(1.5, 0.5);
  Fitter Fit(f_1, vec);

  // Fit.Print("data");

  Fit.fit();
  Fit.GetFitInfo(f_1);

  auto der1 = Fit.Derivate(0, 1);
  auto der2 = Fit.Derivate(0, 2);
  auto der3 = Fit.Derivate(0, 3);

  std::cout << "for x = 0\n"
            << "Image f(0): " << Fit.Image(0)
            << "\n1st Derivative: " << der1.first << "Erro: " << der1.second << "\n2st Derivative: " << der2.first << "Erro: " << der2.second << "\n3st Derivative: " << der3.first << "Erro: " << der3.second << std::endl;

  std::cout << "Integral [0, Pi]: " << Fit.Integrate(0, M_PI) << std::endl;

  // DRAW
  // Fit.DrawFit();
  // Fit.DrawFitErrors();

  return 0;
}

/* RootCalculator::RootCalculator(std::string expressao, double a, double b)
{
  Expressao = expressao;
  TFormula funcao_temp("funcao", Expressao.c_str());
  Funcao = funcao_temp;
  A = a;
  B = b; */

/*
    std::cout << "for x = 0\n"
              << "1st Derivative: " << Fit.Derivate(0, 1).first << "\n2st Derivative: " << Fit.Derivate(0, 2).first << "\n3st Derivative: " << Fit.Derivate(0, 3).first << std::endl;
    std::cout << "Image f(0): " << Fit.Image(0) << std::endl;
    std::cout << "Integral [0,3.14159265]: " << Fit.Integrate(0, 3.14159265) << std::endl;
    */

/*
std::string expression = "[0]+[1]/(x*x-0.028)+[2]/((x*x-0.028)*(x*x-0.028))+[3]*x*x+[4]*x*x*x*x";
// TF1 *f_1 = new TF1("f_diam", expression.c_str(), 0.1, 5);

// MAKE TF1
TF1 *f_1 = new TF1("f", "[0] + [1]*exp(-[2]*x[0])*cos(sqrt(pow([3],2)-pow([2],2))*x + [4])", 0, 10);

// SET PARS
f_1->SetParNames("y0", "y1", "lambda", "w", "p");
f_1->SetParameters(-0.08, 4, 9.3e-4, 1e-2, 8);
*/