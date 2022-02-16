#include "Functor3D.h"

double Functor3D::operator()(double x, double y)
{
    return 0;
};

void Functor3D::Draw(double xi, double xf, double yi, double yf, int num, std::string xtitle, std::string ytitle, std::string ztitle)
{
    TApplication *app = new TApplication("app", nullptr, nullptr);
    c = new TCanvas("canvas", "MyFunction", 0, 0, 1280, 720);

    TRootCanvas *r = (TRootCanvas *)c->GetCanvasImp();
    r->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

    int npx = num;
    int npy = num;
    double x = xi; 
    double y = xi;
    double z;
    int k = 0;
    double dx = (xf - xi) / npx;
    double dy = (yf - yi) / npy;

    TGraph2D *dt = new TGraph2D(npx * npy);
    dt->SetNpy(41);
    dt->SetNpx(40);

    for (int i = 0; i < npx; i++)
    {
        for (int j = 0; j < npy; j++)
        {
            z = operator()(x, y);
            dt->SetPoint(k, x, y, z);
            k++;
            y = y + dy;
        }
        x = x + dx;
        y = yi;
    }
    gStyle->SetPalette(1);
    dt->SetMarkerStyle(20);
    dt->Draw("surf1");
    c->Update();
    app->Run();
};