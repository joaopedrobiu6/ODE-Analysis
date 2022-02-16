#ifndef __OPTICALMAT__
#define __OPTICALMAT__

#include "Material.h"

#include "TF1.h"
#include "TH1F.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TDirectory.h"
#include "TFrame.h"
#include <vector>

class OpticalMat : public Material
{
public:
    OpticalMat(std::string s, Double_t fdens, std::vector<std::pair<float, float>> vec);
    ~OpticalMat() = default;

    void SetRefIndex(std::vector<std::pair<float, float>> points); // pair(wavelength, ref index)
    std::vector<std::pair<float, float>> GetRefIndex();

    void SetFitRefIndex(TF1 *f_); // provide function to be fitted through TF1
    TF1 *GetFitRefIndex();        // return TF1 pointer to fit function

    void DrawRefIndexPoints(); // draw points
    void DrawFitRefIndex();    // draw points and function

    void Print(); // define print for this class

private:
    // method with the fit function
    double FitRefIndex(double *x, double *par);
    // we need to store the refractive index characteric of the material
    std::vector<std::pair<float, float>> refind;
    // we need to store a TF1 pointer to the fit Ref Index function
    TF1 *f;
    Material Mat;
};

#endif