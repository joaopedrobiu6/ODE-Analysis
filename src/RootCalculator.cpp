#include "RootCalculator.h"

RootCalculator::RootCalculator(std::string expressao, double a, double b)
{
  Expressao = expressao;
  TFormula funcao_temp("funcao", Expressao.c_str());
  Funcao = funcao_temp;
  A = a;
  B = b;
};

double RootCalculator::firstDerivative(double x, double h, std::string option)
{
    double Derivative = 0;

    if (option == "fivepoint")
    {
        Derivative = ((Funcao.Eval(x - 2 * h) + 8 * Funcao.Eval(x + h)) - (8 * Funcao.Eval(x - h) + Funcao.Eval(x + 2 * h))) / (12 * h);
    }

    return Derivative;
};

std::vector<double> RootCalculator::Secant_Method()
{
    double a = A;
    std::vector<double> b;
    //QUANTO MENOR O STEP SIZE MAIS ZEROS ENCONTRA
    double step = 0.01;

    double walker = A;
    int numzeros = 0;
    while(walker < B)
    {
        //ENCONTRAR UM INTERVALO QUE CONTÉM O ZERO
        if(Funcao.Eval(walker)*(Funcao.Eval(walker - step)) < 0)
        {
            //ZERO ENTRE WALKER E (WALKER - STEP), IMPLEMENTAR METODO DAS SECANTES
            numzeros++;

            double x = walker - step;
            double epsilon = step;
            double xi;
            double count = 0;

            while(fabs(epsilon) > 0.001 && count < 10)
            {
                xi = x - (Funcao.Eval(x))/firstDerivative(x,0.001,"fivepoint");  
                epsilon = xi - x; 
                count++;
                x = xi;
            }

            if(fabs(xi) < 10e-10)
            {
                b.push_back(0);
            }
            else{
                b.push_back(xi);
            }
        }
        walker += step;
    }

    if(numzeros == 0)
    {
        std::cout << "\nA funcao nao tem zeros" << std::endl; 
        b.push_back(0);
        return b; 
    }
    else
    {
        std::cout << "\nNúmero de Zeros: " << numzeros << std::endl;
        return b;

    }
};

void RootCalculator::Draw(std::vector<double> zeros)
{
    std::cout << "\nZeros Encontrados: " << std::endl;
    for(int i = 0; i < zeros.size();i++)
    {
        std::cout << "Zero[" << i + 1 << "] = " << zeros[i] << std::setprecision(10) << std::endl;
    }

    TApplication app("app", nullptr, nullptr);
    TCanvas *c = new TCanvas("canvas", "Interpolator", 0, 0, 1280, 720);
    c->SetGrid();

    TRootCanvas *r = (TRootCanvas *)c->GetCanvasImp();
    r->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
    
    TF1 *f = new TF1("f",Expressao.c_str(),A,B);

    //GRAFICO PARA FAZER PLOT DOS PONTOS COM LABEL
    auto gr = new TGraph(zeros.size());
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue);

    //TGRAPH COM OS ZEROS DA FUNÇÃO
    for (int i=0;i<zeros.size();i++) {
        gr->SetPoint(i,zeros[i],0);    
    }

    //PARA O EIXO DO TGRAPH SER O MESMO QUE O TF1
    gr->GetYaxis()->SetRangeUser(f->GetMinimum(),f->GetMaximum());
    gr->GetXaxis()->SetRangeUser(A,B);
    gr->Draw("AP");
    f->Draw("same");

    //DESENHAR TEXT COM O VALOR DO ZERO EM CADA ZERO ENCONTRADO
    for(int i = 0; i < zeros.size();i++)
    {
        TText *t = new TText(zeros[i],0,(std::to_string(zeros[i])).c_str()); 
        t->SetTextAlign(11);
        t->SetTextColor(kBlue);
        t->SetTextFont(43);
        t->SetTextSize(20);
        t->SetTextAngle(90);
        t->Draw();

    }

    app.Run();
};