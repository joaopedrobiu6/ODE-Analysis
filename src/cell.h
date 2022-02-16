#ifndef __cell__
#define __cell__
#include <vector>

using namespace std;

//ESTRUTURA CÉLULA
struct cell
{
  int coord_centro [3];          // cm 
  int normal [3] = {0, 0, -1}; // unitary vector
  float vetor_ao_centro [3];
  float area;      // cm^2
  float power;     // W
  float distancia; //m
  float angulo;    //rad
};

#endif
