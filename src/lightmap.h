#ifndef __lightmap__
#define __lightmap__
#include <iostream>
#include "cell.h"
#include <vector>

using namespace std;

class lightmap
{
public:
    // constructors and destructor
    lightmap(const int& ncellx_, const int& ncelly_,vector<vector<cell>> &grid); // number of cells along x and y
    //lightmap(const vector<float> &vx, const vector<float> &vy);
    ~lightmap() = default;

    pair<int, int> GetCellIndex(float x, float y);           // return cell indices
    pair<float, float> GetCellCoo(int index_x, int index_y); // return cell center coo

    double GetCellPower(int index_x, int index_y); // return cell power Watts
    double GetCellPower(float x, float y);         // return cell power Watts

    int GetCellNx(); // get number of cells along x
    int GetCellNy();
    void GetSourceLocation();

    const cell &GetMaxCell();                   // get cell with max power
    std::vector<std::vector<cell>> &GetCells(); // return cells grid

    // (...) other methods you find useful: distance_to_cell(...), ...

private:
    int x,y,z;
    int ncellx;
    int ncelly;
    vector<vector<cell>> GRID;
};

#endif