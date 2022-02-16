#include "Interpolator.h"

Interpolator::Interpolator(int N, double *x1, double *y1) : DataPoints(N, x1, y1)
{
    auto l = [this](double *x, double *p)
    { return this->InterpolateLagrange(x[0]); };
    auto n = [this](double *x, double *p)
    { return this->InterpolateNewton(x[0]); };
    auto sp = [this](double *x, double *p)
    { return this->InterpolateSpline3(x[0]); };
    TF1 *lagrange = new TF1("lagrange", l, 0, 10, 0);
    MI.insert({"lagrange", lagrange});
    TF1 *newton = new TF1("newton", n, 0, 10, 0);
    MI.insert({"newton", newton});
    TF1 *spline3 = new TF1("spline3", sp, 0, 10, 0); 
    MI.insert({"spline3", spline3});
};
Interpolator::Interpolator(const std::vector<std::pair<double, double>> &P1) : DataPoints(P1)
{
    auto l = [this](double *x, double *p)
    { return this->InterpolateLagrange(x[0]); };
    auto n = [this](double *x, double *p)
    { return this->InterpolateNewton(x[0]); };
    auto sp = [this](double *x, double *p)
    { return this->InterpolateSpline3(x[0]); };
    TF1 *lagrange = new TF1("lagrange", l, 0, 10, 0);
    MI.insert({"lagrange", lagrange});
    TF1 *newton = new TF1("newton", n, 0, 10, 0);
    MI.insert({"newton", newton});
    TF1 *spline3 = new TF1("spline3", sp, 0, 10, 0);
    MI.insert({"spline3", spline3});
}; 

Interpolator::Interpolator(const std::vector<std::pair<double, double>> &P1, double min, double max) : DataPoints(P1)
{
    auto l = [this](double *x, double *p)
    { return this->InterpolateLagrange(x[0]); };
    auto n = [this](double *x, double *p)
    { return this->InterpolateNewton(x[0]); };
    auto sp = [this](double *x, double *p)
    { return this->InterpolateSpline3(x[0]); };
    TF1 *lagrange = new TF1("lagrange", l, min, max, 2);
    MI.insert({"lagrange", lagrange});
    TF1 *newton = new TF1("newton", n, min, max, 2);
    MI.insert({"newton", newton});
    TF1 *spline3 = new TF1("spline3", sp, min, max, 2);
    MI.insert({"spline3", spline3});
};

Interpolator::Interpolator(const std::vector<double> &x2, const std::vector<double> &y2, double min, double max) : DataPoints(x2, y2)
{
    auto l = [this](double *x, double *p)
    { return this->InterpolateLagrange(x[0]); };
    auto n = [this](double *x, double *p)
    { return this->InterpolateNewton(x[0]); };
    auto sp = [this](double *x, double *p)
    { return this->InterpolateSpline3(x[0]); };
    TF1 *lagrange = new TF1("lagrange", l, min, max, 0);
    MI.insert({"lagrange", lagrange});
    TF1 *newton = new TF1("newton", n, min, max, 0);
    MI.insert({"newton", newton});
    TF1 *spline3 = new TF1("spline3", sp, min, max, 0);
    MI.insert({"spline3", spline3});
};
Interpolator::~Interpolator()
{
    MI.clear();
    std::cout << "\nCleaning...Done!" << std::endl;
};

TF1 *Interpolator::GetFunction(std::string s)
{
    return MI[s];
};

double Interpolator::InterpolateLagrange(double a)
{
    int N = P.size();
    double yp = 0.;
    double l;

    for (int i = 0; i < N; i++)
    {
        l = 1;
        for (int j = 0; j < N; j++)
        {
            if (i != j)
                l *= (a - P[j].first) / (P[i].first - P[j].first);
        }
        yp += l * P[i].second;
    }
    return yp;
};

void Interpolator::InterpolateNewtonCoefs()
{
    int N = P.size();
    Coefs.resize(N);

    for (int i = 0; i < N; i++)
    {
        Coefs[i] = P[i].second;
    }
    for (int k = 1; k < N; k++)
    {
        for (int i = k; i < N; i++)
        {
            Coefs[i] = (Coefs[i] - Coefs[k - 1]) / (P[i].first - P[k - 1].first);
        }
    }
};

double Interpolator::InterpolateNewton(double r)
{
    int N = P.size();

    double poly = Coefs[N - 1];
    for (int k = 1; k < N; k++)
        poly = Coefs[N - 1 - k] + (r - P[N - 1 - k].first) * poly;

    return poly;
};

double Interpolator::InterpolateSpline3(double a)
{
    int N = P.size();
    int value = 0;
    double poly = 0;

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MK(N - 2, N - 2);
    Eigen::VectorXd b(N - 2);
    Eigen::VectorXd K1(N - 2);
    Eigen::VectorXd K(N);
    MK.setZero();
    b.setZero();
    K1.setZero();
    K.setZero();

    // PREENCHER A MATRIZ MK

    for (int i = 0; i < N - 3; i++)
    {
        MK(i, i) = 2 * (P[i].first - P[i + 2].first);
        MK(i, i + 1) = (P[i + 1].first - P[i + 2].first);
        MK(i + 1, i) = MK(i, i + 1);
    }
    MK(N - 3, N - 3) = 2 * (P[N - 3].first - P[N - 1].first);

    // PREENCHER O VETOR K1

    for (int i = 0; i < N - 2; i++)
    {
        K1[i] = 6 * (((P[i].second - P[i + 1].second) / (P[i].first - P[i + 1].first)) - ((P[i + 1].second - P[i + 2].second) / (P[i + 1].first - P[i + 2].first)));
    }

    // OBTER OS VALORES DO K(1 A N-2)

    EqSolver s3(MK, K1);
    b = s3.GaussSolver(true);

    K[0] = 0;
    K[N - 1] = 0;

    for (int i = 0; i < N - 2; i++)
    {
        K[i + 1] = b[i];
    }

    for (int i = 0; i < N; i++)
    {
        if (P[i].first <= a && P[i + 1].first >= a)
        {
            value = i;
            break;
        }
    }

    if (P[N - 1].first <= a)
    {
        value = N - 2;
    }

    return (K[value] / 6) * ((pow(a - P[value + 1].first, 3)) / (P[value].first - P[value + 1].first) - (a - P[value + 1].first) * (P[value].first - P[value + 1].first)) -
           (K[value + 1] / 6) * ((pow(a - P[value].first, 3)) / (P[value].first - P[value + 1].first) - (a - P[value].first) * (P[value].first - P[value + 1].first)) +
           (P[value].second * (a - P[value + 1].first) - P[value + 1].second * (a - P[value].first)) / (P[value].first - P[value + 1].first);
};

void Interpolator::Draw(std::string s)
{
    TApplication app("app", nullptr, nullptr);
    TCanvas *c = new TCanvas("canvas", "Interpolator", 0, 0, 1280, 720);

    TRootCanvas *r = (TRootCanvas *)c->GetCanvasImp();
    r->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
    DataPoints::Draw();

    if (s == "lagrange")
    {
        MI[s]->SetLineColor(kBlue);
        MI[s]->SetLineWidth(3);
        MI[s]->Draw("SAME");
        // multi_g->Add(MI[s]);
    }
    if (s == "newton")
    {
        this->InterpolateNewtonCoefs();
        MI[s]->SetLineColor(kGreen);
        MI[s]->SetLineWidth(3); 
        MI[s]->Draw("SAME");
        // multi_g->Add(MI[s]);
    }
    if (s == "spline3")
    {
        MI[s]->SetLineColor(kRed);
        MI[s]->SetLineWidth(3);
        MI[s]->Draw("SAME");
        //  multi_g->Add(MI[s]);
    }
    if (s == "all")
    {
        MI["lagrange"]->SetLineColor(kBlue);
        MI["lagrange"]->SetLineWidth(6);
        MI["lagrange"]->SetTitle("Lagrange Method");
        MI["lagrange"]->Draw("SAME");

        this->InterpolateNewtonCoefs();

        MI["newton"]->SetLineColor(kGreen);
        MI["newton"]->SetLineWidth(3);
        MI["newton"]->SetTitle("Newton Method");
        MI["newton"]->Draw("SAME");

        MI["spline3"]->SetLineColor(kRed);
        MI["spline3"]->SetLineWidth(3);
        MI["spline3"]->SetTitle("Spline 3");
        MI["spline3"]->Draw("SAME");

        //c->BuildLegend();
    }
    c->Update();
    c->SaveAs("dist.png");

    app.Run();
};
