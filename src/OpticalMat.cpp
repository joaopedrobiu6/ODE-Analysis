#include "OpticalMat.h"
 
OpticalMat::OpticalMat(std::string s, Double_t fdens, std::vector<std::pair<float, float>> vec)
{
    Material Mat_(s, fdens);
    Mat = Mat_;
    int N = vec.size();
    refind.resize(N);
    for (int i = 0; i < N; i++)
    {
        refind[i].first = vec[i].first;
        refind[i].second = vec[i].second;
    }
    //TF1 *f1 = new TF1("mat", this, &OpticalMat::FitRefIndex, 0, 10, 5);
    //f = f1;
};

/*
OpticalMat::~OpticalMat()
{
    delete f;
};
*/

void OpticalMat::SetRefIndex(std::vector<std::pair<float, float>> points) // pair(wavelength, ref index)
{
    int N = points.size();
    refind.resize(N);

    for (int i = 0; i < N; i++)
    {
        refind[i].first = points[i].first;
        refind[i].second = points[i].second;
    }
};

std::vector<std::pair<float, float>> OpticalMat::GetRefIndex()
{
    return refind;
};

void OpticalMat::SetFitRefIndex(TF1 *f_)
{
    f = f_;
};

TF1 *OpticalMat::GetFitRefIndex()
{
    return f;
};

void OpticalMat::Print()
{
    int N = refind.size();
    Mat.Print();

    // std::cout << " " << std::endl;++++ Material ++++ " << std::endl;Name: " << nome << " " << std::endl;Density: " << densidade << std::endl;
    std::cout << "(Wavelenght, Ref Index)"
              << std::endl;
    for (int i = 0; i < N; i++)
    {
        std::cout << "(" << refind[i].first << ", " << refind[i].second << ")" << std::endl;
    }
};

double OpticalMat::FitRefIndex(double *x, double *par)
{
    double ind;
    ind = par[0] + (par[1] / (x[0] - 0.028)) + (par[2] / pow(x[0] * x[0] - 0.028, 2)) + par[3] * x[0] * x[0] + par[4] * x[0] * x[0] * x[0] * x[0];
    return ind;
};

void OpticalMat::DrawRefIndexPoints()
{
    TCanvas *c = new TCanvas("c", "Graph Refraction Index", 1280, 720);
    c->SetGrid();
    c->GetFrame()->SetFillColor(21);
    c->GetFrame()->SetBorderSize(12);
    int n = refind.size();
    TGraph *gr = new TGraph();

    for (int i = 0; i < n; i++)
    {
        gr->SetPoint(i, refind[i].first, refind[i].second);
    }

    gr->SetTitle("Refraction Index");
    gr->GetXaxis()->SetTitle("Wavelenght [nm]");
    gr->GetYaxis()->SetTitle("Ref Index");
    gr->SetMarkerColor(kBlue);
    gr->SetMarkerStyle(21);
    gr->Draw("AP");
    c->SaveAs("points.pdf");
};

void OpticalMat::DrawFitRefIndex()
{
    TCanvas *c1 = new TCanvas("c1", "Graph Ref_Index", 1280, 720);
    c1->SetGrid();
    int n = refind.size();
    TGraph *gr1 = new TGraph();

    for (int i = 0; i < n; i++)
    {
        gr1->SetPoint(i, refind[i].first, refind[i].second);
    }

    gr1->Fit(f);

    gr1->SetTitle("Refraction Index");
    gr1->GetXaxis()->SetTitle("Wavelenght [nm]");
    gr1->GetYaxis()->SetTitle("Ref Index");

    gr1->SetLineColor(69);
    gr1->SetLineWidth(4);
    gr1->SetMarkerColor(kBlue);
    gr1->SetMarkerStyle(21);
    gr1->Draw("AP");
    c1->SaveAs("fit.pdf");
};