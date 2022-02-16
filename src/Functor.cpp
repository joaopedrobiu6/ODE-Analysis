#include "Functor.h"

double Functor::operator()(double x)
{
    return 0;
}; 

void Functor::Draw(double xi, double xf, int num, std::string xtitle, std::string ytitle, std::string option)
{
    double_t x_values[num]; // create vector of n values
    double_t y_values[num];
    double x_inicial = xi, step = (xf - xi) / num;

    // Encher vector
    for (int i = x_inicial; i < num; i++)
    {
        x_values[i] = x_inicial + i * step;
        y_values[i] = operator()(x_values[i]); 
    }

    TApplication *app = new TApplication("app", nullptr, nullptr);
    c = new TCanvas("canvas", "MyFunction", 0, 0, 1280, 720);

    TRootCanvas *r = (TRootCanvas *)c->GetCanvasImp();
    r->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

    TGraph *gr = new TGraph(num, x_values, y_values);
    gROOT->SetStyle("Plain"); 
    gr->SetLineColor(69);
    gr->SetLineWidth(4);
    gr->SetName("Graph");
    gr->GetXaxis()->SetTitle(xtitle.c_str());
    gr->GetYaxis()->SetTitle(ytitle.c_str());
    gr->SetTitle("Function");
    gr->SetFillColor(97);

    if (option == "Integral")
        gr->Draw("ABC");
    else
        gr->Draw("AC");

    c->Update();
    app->Run();
};
