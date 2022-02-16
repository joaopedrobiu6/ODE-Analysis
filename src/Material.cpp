#include "Material.h"

void Material::Print()
{
    std::string nome = GetName();
    Double_t densidade = GetDensity();

    std::cout << "\n###  Material  ###\nName: " << nome << "\nDensity: " << densidade << std::endl;
};
