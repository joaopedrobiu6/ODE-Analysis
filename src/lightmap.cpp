 #include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include <cmath>
#include <math.h>
#include "lightmap.h"
#include "cell.h"

using namespace std;

// Construtor
lightmap::lightmap(const int& ncellx_, const int& ncelly_,vector<vector<cell>> &grid)
{
    ncellx = ncellx_;
    ncelly = ncelly_;
    GRID = grid;

}; // number of cells along x and y

//lightmap::lightmap(const vector<float> &vx, const vector<float> &vy);

void lightmap::GetSourceLocation()
{
    //INPUT DAS COORDENADAS DA FONTE DE LUZ
    cout << "Inserir Coordenadas da fonte de luz" << endl;
    cout << "x:";
    cin >> x;
    cout << "y:";
    cin >> y;
    cout << "z:";
    cin >> z;

}

double lightmap::GetCellPower(int index_x, int index_y)
{
    //POTÊNCIA
    const int fluxo_radiante = 100; //watts

    //COORDENADAS DO CENTRO DAS GRIDS
    GRID[index_x][index_y].coord_centro[0] = index_x - 150; //x
    GRID[index_x][index_y].coord_centro[1] = index_y - 100; //y
    GRID[index_x][index_y].coord_centro[2] = 100;           //z

    //ÁREA DAS CÉLULAS
    GRID[index_x][index_y].area = (300 / ncellx) * (200 / ncelly); //cm^2

    //COORDENADAS DO VETOR CENTRO_CÉLULA - FONTE
    GRID[index_x][index_y].vetor_ao_centro[0] = GRID[index_x][index_y].coord_centro[0] - x;
    GRID[index_x][index_y].vetor_ao_centro[1] = GRID[index_x][index_y].coord_centro[1] - y;
    GRID[index_x][index_y].vetor_ao_centro[2] = GRID[index_x][index_y].coord_centro[2] - z;

    //NORMA DO VETOR DO CENTRO DA CÉLULA À FONTE - RAIO
    GRID[index_x][index_y].distancia = sqrt((GRID[index_x][index_y].coord_centro[0] - x) * (GRID[index_x][index_y].coord_centro[0] - x) + (GRID[index_x][index_y].coord_centro[1] - y) * (GRID[index_x][index_y].coord_centro[1] - y) + (GRID[index_x][index_y].coord_centro[2] - z) * (GRID[index_x][index_y].coord_centro[2] - z));

    double distanc = sqrt((GRID[index_x][index_y].vetor_ao_centro[0] * GRID[index_x][index_y].vetor_ao_centro[0]) + (GRID[index_x][index_y].vetor_ao_centro[1] * GRID[index_x][index_y].vetor_ao_centro[1]) + (GRID[index_x][index_y].vetor_ao_centro[2] * GRID[index_x][index_y].vetor_ao_centro[2]));

    //ÂNGULO ENTRE O VETOR NORMAL DA CÉLULA E O VETOR À FONTE DE LUZ
    GRID[index_x][index_y].angulo = acos((GRID[index_x][index_y].vetor_ao_centro[2]) / distanc);

    //CÁLCULO DA POTÊNCIA
    double power = abs(fluxo_radiante * (cos(GRID[index_x][index_y].angulo))) / (4 * M_PI * (pow(distanc, 2) / 10000)); //watt/cm^2;

    return power;
};

int lightmap::GetCellNx()
{
    cout << "\nNumber of cells - xx:" << ncellx << endl;
    return ncellx;
};
int lightmap::GetCellNy()
{
    cout << "\nNumber of cells - yy:" << ncelly << endl;
    return ncelly;
};
pair<int, int> GetCellIndex(float x, float y){ //ACABAR
pair<int, int> v = {1, 1};
return v;
};
pair<float, float> GetCellCoo(int index_x, int index_y){ //ACABAR
pair<float, float> v = {1, 1};

return v;
}