#include "FourierDT.h"

std::vector<std::pair<double, double>> FourierDT::Fourier()
{
    std::cout << "\n"
              << __FUNCTION__ << std::endl;

    std::vector<std::pair<double, double>> X;
    double N = resultado.size();
    double sum, sum_im;
    double f_k = 0;
    int k = 0;
    while (f_k < 5)
    {
        sum = 0;
        sum_im = 0;
        f_k = (double)k / ((double)N * step); // tc
        for (int n = 0; n < N; n++)
        {
            sum += resultado[n].X()[0] * (cos(2 * M_PI / N * n * k));    // Parte Real
            sum_im += resultado[n].X()[0] * (sin(2 * M_PI / N * n * k)); // Parte Imaginária
        }
        X.push_back(std::make_pair(f_k, sqrt(sum * sum + sum_im * sum_im) / N));
        k++;
    }

    return X;
};

void FourierDT::Fourier_Draw_Harmonics(const char *s, std::vector<std::pair<double, double>> Xk)
{
    std::cout << "\n"
              << __FUNCTION__ << std::endl;

    int N = Xk.size();

    TGraph *gr = new TGraph(N);

    for (int k = 0; k < N; k++)
    {
        gr->SetPoint(k, Xk[k].first, Xk[k].second);
    }

    TCanvas *c = new TCanvas("canvas", "pendulum_parametric", 0, 0, 3200, 1800);
    c->SetTickx();
    c->SetTicky();
    c->SetGridx();
    c->SetGridy();

    // DESENHAR AS 5 FREQUÊNCIAS MAIS IMPORTANTES
    std::vector<double> max;
    std::vector<double> freq;
    std::cout << "\nFrequências Fundamentais:" << std::endl;
    for (int i = 1; i <= N - 1; i++)
    {
        if (Xk[i].second > Xk[i - 1].second && Xk[i].second > Xk[i + 1].second)
        {
            max.push_back(Xk[i].second);
            freq.push_back(Xk[i].first);
            std::cout << "Máximo em: " << Xk[i].first << "Hz, Amplitude: " << Xk[i].second << std::endl;
        }
    }
    std::cout << std::endl;

    gr->SetTitle("|X(f_{k})| vs. f_{k}");
    gr->GetXaxis()->SetTitle("f_{k} [Hz]");
    gr->GetYaxis()->SetTitle("|X(f_{k})|");
    gr->SetFillColor(50);

    gStyle->SetOptStat(0);

    // h->Draw("HIST");
    gr->Draw("AB + same");

    for (int i = 0; i < freq.size(); i++)
    {
        std::string max_str = std::to_string(freq[i]) + " Hz";
        TText *t = new TText(freq[i], max[i], max_str.c_str());
        t->SetTextAlign(11);
        t->SetTextColor(1);
        t->SetTextFont(43);
        t->SetTextSize(55);
        t->SetTextAngle(90);
        t->Draw();
    }

    c->Update();
    c->SaveAs(s);

    delete c;
    delete gr;
};

void FourierDT::Fourier_Draw_Max(const char *s, std::vector<std::pair<double, double>> Xk, int qual_deles)
{
    std::cout << "\n"
              << __FUNCTION__ << std::endl;

    int N = Xk.size();

    TGraph *gr = new TGraph(N);

    for (int k = 0; k < N; k++)
    {
        gr->SetPoint(k, Xk[k].first, Xk[k].second);
    }

    TCanvas *c = new TCanvas("canvas", "pendulum_parametric", 0, 0, 3200, 1800);
    c->SetTickx();
    c->SetTicky();
    c->SetGridx();
    c->SetGridy();

    // DESENHAR AS 5 FREQUÊNCIAS MAIS ALTAS
    std::vector<double> max;
    std::vector<double> freq;
    std::cout << "\nMáximos obtidos:" << std::endl;
    for (int k = 0; k < 5; k++)
    {
        double max_temp = -1;
        int max_position = 0;
        for (int i = 0; i < N; i++)
        {
            if (Xk[i].second > max_temp)
            {
                max_temp = Xk[i].second;
                max_position = i;
            }
        }
        max.push_back(Xk[max_position].second);
        freq.push_back(Xk[max_position].first);
        std::cout << "Máximo em: " << Xk[max_position].first << "Hz, Amplitude: " << Xk[max_position].second << std::endl;
        Xk[max_position].second = 0;
    }
    std::cout << std::endl;

    gr->SetTitle("|X(f_{k})| vs. f_{k}");
    gr->GetXaxis()->SetTitle("f_{k} [Hz]");
    gr->GetYaxis()->SetTitle("|X(f_{k})|");
    gr->SetFillColor(50);

    gStyle->SetOptStat(0);

    // h->Draw("HIST");
    gr->Draw("AB + same");

    int a; // Para não sobrepor valores no grafico
    if (qual_deles != 0)
    {
        a = freq.size();
    }
    else
        a = 4;
    for (int i = 0; i < a; i++)
    {
        std::string max_str = std::to_string(freq[i]) + " Hz";
        TText *t = new TText(freq[i], max[i], max_str.c_str());
        t->SetTextAlign(11);
        t->SetTextColor(1);
        t->SetTextFont(43);
        t->SetTextSize(55);
        t->SetTextAngle(qual_deles);
        t->Draw();
    }

    c->Update();
    c->SaveAs(s);

    delete c;
    delete gr;
};

void FourierDT::Comparacao_analitica(const char *filename, const char *xaxis, const char *yaxis, int color, std::vector<std::pair<double, double>> Freq)
{
    std::cout << "\n"
              << __FUNCTION__ << " - Imprimindo o gráfico da soluções da ODE..." << std::endl;
    // FAZER GRÁFICO DA EDO
    TCanvas *c = new TCanvas("canvas", "pendulum_parametric", 0, 0, 3200, 1800);
    c->SetTickx();
    c->SetTicky();
    c->SetGridx();
    c->SetGridy();

    int n = resultado.size();

    double *tempo = new double[n];
    double *posicao = new double[n];

    for (int i = 0; i < n; i++)
    {
        tempo[i] = resultado[i].T();
        posicao[i] = resultado[i].X()[0];
    }

    TGraph *gr = new TGraph(n, tempo, posicao);
    gStyle->SetTitleFontSize(0.05);
    gr->SetLineColor(50);
    gr->SetLineWidth(10);

    gr->SetTitle("Aproximacao da Solucao Numerica a uma Solucao Explicita");
    gr->GetXaxis()->SetTitle(xaxis);
    gr->GetYaxis()->SetTitle(yaxis);

    // Encontrar as harmonicas
    std::vector<double> max;
    std::vector<double> freq;
    for (int i = 1; i <= Freq.size() - 1; i++)
    {
        if (Freq[i].second > Freq[i - 1].second && Freq[i].second > Freq[i + 1].second)
        {
            max.push_back(Freq[i].second);
            freq.push_back(Freq[i].first);
        }
    }

    TF1 *f = new TF1("fa2", "2*([0]*cos(2*TMath::Pi()*[1]*(x - 6)) + [2]*cos(2*TMath::Pi()*[3]*(x - 6)) +[4]*cos(2*TMath::Pi()*[5]*(x - 6)) +[6]*cos(2*TMath::Pi()*[7]*(x - 6))) ", 0, 100);
    f->SetParameter(0, max[0]);
    f->SetParameter(1, freq[0]);
    f->SetParameter(2, max[1]);
    f->SetParameter(3, freq[1]);
    f->SetParameter(4, max[2]);
    f->SetParameter(5, freq[2]);
    f->SetParameter(6, max[3]);
    f->SetParameter(7, freq[3]);

    f->SetLineColor(52);
    f->SetLineWidth(10);
    gr->Draw("AC");
    f->Draw("AC + same");

    auto legenda = new TLegend(0.7, 0.25, 0.9, 0.1);
    legenda->SetFillStyle(0);

    legenda->AddEntry(gr, "Runge-Kutta 4", "L");
    legenda->AddEntry(f, "Aproximacao Explicita", "L");
    legenda->Draw();

    c->SaveAs(filename);

    delete[] tempo;
    delete[] posicao;

    delete c;
};
