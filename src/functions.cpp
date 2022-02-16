#include "functions.h"

functions::functions(double g_, double L_)
{
    g = g_;
    L = L_;
};

double functions::get_g()
{
    return g;
};

double functions::get_L()
{
    return L;
};

void functions::ODE_Draw(const char *filename, const char *xaxis, const char *yaxis, int color, std::vector<ODEpoint> resultado)
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
    gr->SetLineColor(color);
    gr->SetLineWidth(10);

    gr->SetTitle("Pendulo");
    gr->GetXaxis()->SetTitle(xaxis);
    gr->GetYaxis()->SetTitle(yaxis);
    gr->Draw("AC");

    c->SaveAs(filename);

    delete[] tempo;
    delete[] posicao;

    delete c;
};

std::vector<double> functions::Moving_Average(std::vector<ODEpoint> vetor, double time_step, double tw)
{
    double sum_sq;
    std::vector<double> amp_mean;
    amp_mean.resize(vetor.size() - tw);
    int loop;

    for (int i = 0; i < vetor.size() - tw; i++)
    {
        amp_mean[i] = 0;
    }

    loop = 0;
    while (loop < tw)
    {
        sum_sq += vetor[loop].X()[0] * vetor[loop].X()[0];
        loop++;
    }
    double media = sqrt(sum_sq / tw);
    amp_mean[0] = media;

    for (int i = 1; i < vetor.size() - tw; i++)
    {
        sum_sq += vetor[i + tw].X()[0] * vetor[i + tw].X()[0];
        sum_sq -= vetor[i - 1].X()[0] * vetor[i - 1].X()[0];
        media = sqrt(sum_sq / tw);
        amp_mean[i] = media;
    }

    return amp_mean;
};

void functions::Moving_Average_Draw(const char *filename, const char *xaxis, const char *yaxis, int color, std::vector<ODEpoint> resultado, std::vector<double> media, double epsilon, double time_step, double tw)
{
    std::cout << "\n"
              << __FUNCTION__ << " - Imprimindo Moving Average..." << std::endl;
    // FAZER GRÁFICO DA EDO
    TCanvas *c = new TCanvas("canvas", "pendulum_parametric", 0, 0, 3200, 1800);
    c->SetTickx();
    c->SetTicky();
    c->SetGridx();
    c->SetGridy();

    int n = resultado.size();

    double *tempo = new double[n];
    double *posicao = new double[n];
    double *velocidade = new double[n];

    for (int i = 0; i < n; i++)
    {
        tempo[i] = resultado[i].T();
        posicao[i] = resultado[i].X()[0];
        velocidade[i] = resultado[i].X()[1];
    }

    TGraph *gr = new TGraph(n, tempo, posicao);
    gStyle->SetTitleFontSize(0.05);
    gr->SetLineColor(50);
    gr->SetLineWidth(10);

    gr->SetTitle("Grafico Media Deslizante");
    gr->GetXaxis()->SetTitle(xaxis);
    gr->GetYaxis()->SetTitle(yaxis);
    gr->Draw("AC");

    double *media_arr = &media[0];

    double *tempo_novo = new double[n - 5000]; // Array com novo tamanho para descartar as medias que não foram calculadas nos extremos

    for (int i = 0; i < n - tw; i++)
    {
        tempo_novo[i] = tempo[(int)(i + tw / 2)];
    }

    TGraph *gr2 = new TGraph(n - tw, tempo_novo, media_arr);
    gr2->SetLineColor(62);
    gr2->SetLineWidth(10);
    gr2->Draw("same");

    auto legenda = new TLegend(0.7, 0.25, 0.9, 0.1);
    legenda->SetFillStyle(0);

    std::string s_epsilon;
    if (epsilon == 0.36)
        s_epsilon = "0.36";
    else
        s_epsilon = "0.2";

    std::string legenda_texto = "#delta = 0.4 #varepsilon = " + s_epsilon; // Criar string para a legenda com o valor da particula
    legenda->AddEntry(gr, legenda_texto.c_str(), "L");
    legenda->AddEntry(gr2, "Media Deslizante", "L");
    legenda->Draw();
    c->SaveAs(filename);

    delete[] tempo;
    delete[] posicao;
    delete[] velocidade;
    delete[] tempo_novo;

    delete c;
};

void functions::Phase_Draw(const char *title, const char *s, int color, std::vector<ODEpoint> resultado)
{
    std::cout << "\n"
              << __FUNCTION__ << " - Imprimindo o Retrato de Fase..." << std::endl;
    // FAZER GRÁFICO DA EDO
    TCanvas *c = new TCanvas("canvas", "pendulum_parametric", 0, 0, 6100, 3600);
    c->SetTickx();
    c->SetTicky();
    c->SetGridx();
    c->SetGridy();

    int n = resultado.size();

    double posicao[n];
    double velocidade[n];

    for (int i = 0; i < n; i++)
    {
        posicao[i] = resultado[i].X()[0];
        velocidade[i] = resultado[i].X()[1];
    }

    TGraph *gr = new TGraph(n, posicao, velocidade);
    gStyle->SetTitleFontSize(0.05);
    gr->SetLineColor(color);
    gr->SetLineWidth(10);

    gr->SetTitle(title);
    gr->GetYaxis()->SetTitle("#dot{#hat{#theta}}(#hat{t})");
    gr->GetXaxis()->SetTitle("#hat{#theta}(#hat{t})");

    gr->Draw("AC");

    c->SaveAs(s);

    delete c;
};

int functions::Is_Stable(double *media, int n)
{
    int flag_estab = 1;
    for (int i = 0; i < n; i++)
    {
        if (media[i] > M_PI)
        {
            flag_estab = 0;
            break;
        }
    }
    return flag_estab;
};

void functions::Draw_Stability(const char *filename, int numero, double step)
{
    std::cout << "\n"
              << __FUNCTION__ << " - Imprimindo o Diagrama de Estabilidade..." << std::endl;

    int n;
    int count = 0;
    double step_delta = (1.2 - 0.02) / numero;
    double step_epsilon = (0.5 - 0) / numero;

    auto h = new TH2D("h", "Stability Regions", numero + 1, 0.02, 1.2, numero + 1, 0, 0.5);

    double est = 1;

    for (double delta = 0.02; delta <= 1.2; delta = delta + step_delta)
    {
        for (double epsilon = 0; epsilon <= 0.5 + step_epsilon; epsilon = epsilon + step_epsilon)
        {
            pendulum_parametric pendulum(L, {0.1, 0});
            pendulum.SetFunction(0, [](ODEpoint p)
                                 { return p.X()[1]; });

            pendulum.SetFunction(1, [&](ODEpoint p)
                                 { return -((delta + epsilon * std::cos(p.T())) * std::sin((p.X())[0])); });

            std::vector<ODEpoint> resultado = pendulum.LeapFrogImprovedSolver(100, step);
            std::vector<double> mov_avg = Moving_Average(resultado, 5e-3, 5 / 5e-3);

            n = mov_avg.size();

            double *media = &mov_avg[0];

            int value = Is_Stable(media, n);
            if (value == 0)
                est = 0.1;
            else
                est = 1;

            h->Fill(delta, epsilon, est);
        }
        count++;
        if (count % 5 == 0)
            std::cout << 100 * delta / 1.2 << "% " << std::endl;
    }
    std::cout << std::endl;

    TCanvas *c = new TCanvas("canvas", "pendulum_parametric", 0, 0, 6400, 3600);

    TStyle *mcStyle = new TStyle("mcStyle", "FC's Root Styles");
    mcStyle->SetPalette(kGreyYellow);

    h->SetTitle("Diagrama de Estabilidade para (#delta, #varepsilon)");
    h->GetXaxis()->SetTitle("#delta");
    h->GetYaxis()->SetTitle("#varepsilon");

    gStyle->SetOptStat(0);

    h->SetContour(1000);
    h->SetStats(0);
    h->SetLineColor(1);
    h->SetBarWidth(2);
    h->Draw("COL");

    // TGraph para as cores
    TGraph *gr1 = new TGraph();
    TGraph *gr2 = new TGraph();
    gr1->SetLineColor(92);
    gr2->SetLineColor(112);
    gr1->SetLineWidth(20);
    gr2->SetLineWidth(20);
    // gr1->Draw("same");
    // gr2->Draw("same");

    auto legenda = new TLegend(0.7, 0.25, 0.9, 0.1);
    legenda->SetFillColor(0);
    legenda->AddEntry(gr1, "Solucoes Estaveis", "LP");
    legenda->AddEntry(gr2, "Solucoes Instaveis", "LP");
    legenda->Draw();

    c->Update();
    c->SaveAs(filename);

    delete c;
    delete h;
};

void functions::Draw_Stability_Zoom(const char *filename, int numero, double step)
{
    std::cout << "\n"
              << __FUNCTION__ << " - Imprimindo o Diagrama de Estabilidade Aumentado (Zoom In)..." << std::endl;

    int count = 0;
    int n;
    double step_delta = (0.35 - 0.2) / numero;
    double step_epsilon = (0.25 - 0.1) / numero;

    auto h = new TH2D("h", "Stability Regions", numero + 1, 0.2, 0.35, numero + 1, 0.1, 0.25);

    double est = 1;

    for (double delta = 0.2; delta <= 0.35; delta = delta + step_delta)
    {
        for (double epsilon = 0.1; epsilon <= 0.25 + step_epsilon; epsilon = epsilon + step_epsilon)
        {
            pendulum_parametric pendulum(L, {0.1, 0});
            pendulum.SetFunction(0, [](ODEpoint p)
                                 { return p.X()[1]; });

            pendulum.SetFunction(1, [&](ODEpoint p)
                                 { return -((delta + epsilon * std::cos(p.T())) * std::sin((p.X())[0])); });

            std::vector<ODEpoint> resultado = pendulum.LeapFrogImprovedSolver(100, step);
            std::vector<double> mov_avg = Moving_Average(resultado, step, 5 / step);

            n = mov_avg.size();

            double *media = &mov_avg[0];

            int value = Is_Stable(media, n);
            if (value == 0)
                est = 0.1;
            else
                est = 1;

            h->Fill(delta, epsilon, est);
        }
        count++;
        if (count % 5 == 0)
            std::cout << 100 * delta / 0.35 << "% " << std::endl;
    }

    TCanvas *c = new TCanvas("canvas", "pendulum_parametric", 0, 0, 6400, 3600);

    TStyle *mcStyle = new TStyle("mcStyle", "FC's Root Styles");
    mcStyle->SetPalette(kGreyYellow);

    h->SetTitle("Diagrama de Estabilidade para (#delta, #varepsilon)");
    h->GetXaxis()->SetTitle("#delta");
    h->GetYaxis()->SetTitle("#varepsilon");

    gStyle->SetOptStat(0);

    h->SetContour(1000);
    h->SetStats(0);
    h->SetLineColor(1);
    h->SetBarWidth(2);
    h->Draw("COL");

    // TGraph para as cores
    TGraph *gr1 = new TGraph();
    TGraph *gr2 = new TGraph();
    gr1->SetLineColor(92);
    gr2->SetLineColor(112);
    gr1->SetLineWidth(20);
    gr2->SetLineWidth(20);
    // gr1->Draw("same");
    // gr2->Draw("same");

    auto legenda = new TLegend(0.7, 0.25, 0.9, 0.1);
    legenda->SetFillColor(0);
    legenda->AddEntry(gr1, "Solucoes Estaveis", "LP");
    legenda->AddEntry(gr2, "Solucoes Instaveis", "LP");
    legenda->Draw();

    c->Update();
    c->SaveAs(filename);

    delete c;
};

std::vector<std::pair<double, double>> functions::Fourier(std::vector<ODEpoint> resultado, double step)
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

void functions::Fourier_Draw_Harmonics(const char *s, std::vector<std::pair<double, double>> Xk)
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

void functions::Fourier_Draw_Max(const char *s, std::vector<std::pair<double, double>> Xk, int qual_deles)
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

void functions::Comparacao_analitica(const char *filename, const char *xaxis, const char *yaxis, int color, std::vector<ODEpoint> resultado, std::vector<std::pair<double, double>> Freq)
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
