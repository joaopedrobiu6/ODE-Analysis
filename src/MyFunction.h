#ifndef __MYFUNCTION__
#define __MYFUNCTION__

#include "Functor.h"

class MyFunction : public Functor
{
public:
    MyFunction(std::string s = "Functor") : Functor(s) {;}
    ~MyFunction() = default;
    double operator()(double x); 
};
#endif 