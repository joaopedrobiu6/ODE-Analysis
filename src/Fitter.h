#ifndef __FITTER__
#define __FITTER__

#include "TAxis.h"
#include "TROOT.h"
#include "TF1.h"
#include "TH2F.h"
#include "TApplication.h"
#include "TRootCanvas.h"
#include <TROOT.h>
#include <TStyle.h>
#include "TH1F.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TDirectory.h"
#include "TFrame.h"
#include "TGraphErrors.h"

#include <vector>
#include <iostream>
#include <sstream> // std::stringstream
#include <fstream> // std::fstream
#include <iomanip>

#include "DataReader.h"

class Fitter  
{
public:
    // Construtores e Destrutor
    Fitter();
    Fitter(std::vector<std::vector<float>> vec) : data(vec) { ; }
    Fitter(TF1 *f_, std::vector<std::vector<float>> vec) : f(f_), data(vec) { ; }
    Fitter(TF1 *f_) : f(f_) { ; }
    ~Fitter() = default;

    //DATA INPUTS AND GETTERS
    void InputData(std::vector<std::vector<float>> vec); // inserir dados de vec no data (private)
    std::vector<std::vector<float>> GetData();           // retorna data (private) para fora
    void SetData();                                      // ler dados de um ficheiro txt (podiamos mudar para csv...R style)
    void SetDataFromDataReader(std::string);             // ler dados usando DataReader.h

    //FIT
    void fit();
    void SetFit(TF1 *f_);                                // criar TF1 f a partir de f_ (private)
    TF1 *GetFit();                                       // retorna f (private) para fora
    void GetFitInfo(TF1 *f_);                            // imprime os dados do fit, nome das variáveis livres e os seus valores, chi squared, etc

    //DRAW
    void DrawPoints();                                   // desenha os pontos só
    void DrawFitErrors();                                // faz o fit e desenha, com as considerações dos erros ex ey
    void DrawFit();                                      // faz o fit e desenha, sem os erros
    
    //EXTRAS
    std::pair<double, double> Derivate(double x, double order);
    double Image(double x);
    double Integrate(double min, double max);

    //DUMP
    void Print(std::string s = "data");                  // imprime data

private:
    std::vector<std::vector<float>> data;
    TF1 *f;
};

#endif