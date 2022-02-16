#include "MyFunction3D.h"

double MyFunction3D::operator()(double x, double y)
{
    return x/sqrt((x*x+y*y));
};