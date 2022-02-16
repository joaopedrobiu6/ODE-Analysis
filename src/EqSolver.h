#ifndef __EqSolver__
#define __EqSolver__
#include <Eigen/Core>
#include <iostream>
#include "FCmatrixAlgo.h"

using namespace std;

class EqSolver
{
public:
    // constructors and destructor
    EqSolver();
    EqSolver(
        const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &, // matrix coeffs
        const Eigen::VectorXd &);
    ~EqSolver() = default;

    // output (optional)
    friend ostream &operator<<(ostream &, const EqSolver &);

    // solvers
    const Eigen::VectorXd &GaussSolver(bool pivot);
    const Eigen::VectorXd &LUSolver(bool pivot);
    void IterativeJacobiSolver(Eigen::VectorXd &,
                               const int &itmax,
                               double tol);
    void IterativeJacobiSolverRelax(Eigen::VectorXd &,
                                    const int &,
                                    double,
                                    int,
                                    int);
    void IterativeGaussSeidelSolver(
        Eigen::VectorXd &,
        int &itmax,
        double tol);

private:
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> M; // coefficients matrix
    Eigen::VectorXd b;                                       // constants vector
    Eigen::VectorXd sol1;
};

#endif
