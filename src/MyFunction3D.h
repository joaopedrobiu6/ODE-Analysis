#ifndef __MYFUNCTION3D__
#define __MYFUNCTION3D__

#include "Functor3D.h"

class MyFunction3D : public Functor3D
{
public:
    MyFunction3D(std::string s = "Functor3D") : Functor3D(s) {;}
    ~MyFunction3D() = default;
    double operator()(double x, double y);
};
#endif 