#include <iostream>
#include <vector>
#include "DataPoints.h"

DataPoints::DataPoints(const std::vector<std::pair<double, double>> &vetor)
{
    P = vetor;
    graph = new TGraph(P.size());
};

DataPoints::DataPoints(int n, double *x1, double *y1)
{
    int N = n;
    P.resize(N);
    for (int i = 0; i < N; ++i)
    {
        P[i].first = x1[i];
        P[i].second = y1[i];
    }
    graph = new TGraph(P.size());
};

DataPoints::DataPoints(const std::vector<double> &vecx_, const std::vector<double> &vecy_)
{
    std::pair<double, double> par;
    for (int i = 0; i < vecx_.size(); i++)
    {
        par = std::make_pair(vecx_[i], vecy_[i]);
        P.push_back(par);
    }
    graph = new TGraph(P.size());
};

const std::vector<std::pair<double, double>> &DataPoints::GetPoints()
{
    return P;
};

void DataPoints::GetGraph(TGraph &graph)
{
    for (int i = 0; i < P.size(); ++i)
    {
        graph.SetPoint(i, P[i].first, P[i].second);
    }
};

TGraph &DataPoints::GetGraphPointer()
{
    return *graph;
};

void DataPoints::Draw()
{
    GetGraph(*graph);

    graph->SetTitle();
    graph->Draw("AP");
    gROOT->SetStyle("Plain");
    graph->SetName("Graph");
    graph->SetMarkerStyle(kFullCircle);
    graph->SetMarkerSize(1.5);
    graph->SetMarkerColor(kBlack);
    graph->GetXaxis()->SetTitle("x");
    graph->GetYaxis()->SetTitle("y");
    //graph->GetYaxis()->SetRangeUser(-7, -3.5);

    graph->SetTitle("Interpolator");
};

std::ostream &operator<<(std::ostream &s, const DataPoints &A)
{
    auto R = A;
    s << "Número de Pontos: " << R.GetPoints().size() << std::endl;

    for (int i = 0; i < R.GetPoints().size(); ++i)
    {
        s << "(" << R.GetPoints()[i].first << "," << R.GetPoints()[i].second << ")" << std::endl;
    }
    return s;
};