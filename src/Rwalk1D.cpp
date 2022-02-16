#include "Rwalk1D.h"

Rwalk1D::Rwalk1D(int N_, double x_, double pL_, double pR_, double dt_, double dx_)
{
    N = N_;
    x0 = x_;
    pL = pL_;
    pR = pR_;
    dt = dt_;
    dx = dx_;

};

Rwalk1D::~Rwalk1D()
{
    mT.clear();
};

void Rwalk1D::Run(int nsteps)
{
    srand(time(NULL));

    mT.clear();

    std::vector<double> temp;
    double pos_atual;
    double prob;

    for(int i = 0; i < N;i++)
    {
        pos_atual = x0;
        for(int k = 0;k <= nsteps;k++)
        {
            prob = std::rand() % 100;
            if(prob <  pR*100)//Se prob inferior à pR mover para a direita
            {
                pos_atual--;
                temp.push_back(pos_atual);
            } 
            else if(prob >= pR*100 && prob <= (pR*100 + pL*100))//esquerda
            {
                pos_atual++;
                temp.push_back(pos_atual);
            }
            else {
                temp.push_back(pos_atual);
            }
            //std::cout << k << " " << "prob: " << prob << " " << "pos: " << pos_atual << std::endl;
        }
        mT[i].swap(temp);
    }

};

// Getters
const std::vector<double>& Rwalk1D::GetTrajectory(int n)
{
    return mT[n];
};

double Rwalk1D::GetTimeStep(){
    return dt;
};

double Rwalk1D::GetSpaceStep(){
    return dx;
};TMultiGraph *mg = new TMultiGraph();