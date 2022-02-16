#ifndef __FCmatrixAlgo__
#define __FCmatrixAlgo__
#include <Eigen/Core>
#include <iostream>

using namespace std;

class FCmatrixAlgo
{
public:
    FCmatrixAlgo() = default; // compiler do it
    ~FCmatrixAlgo() = default;
    /*
    Implements Gauss elimination
    */
    static void GaussElimination(
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &, // matrix coeffs
        Eigen::VectorXd & // vector of constants
    ); // no pivoting
    static void GaussEliminationPivot(
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &, // matrix coeff
        Eigen::VectorXd &,// vector of constants
        Eigen::VectorXd & // row order indexing
    ); // make pivoting
    /*
    Implements LU decomposition (Doolitle)
    */
    static void LUdecomposition(
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &, // matrix coeff
        Eigen::VectorXd &
    );// no pivoting
    static void LUdecompositionPivot(
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &, // matrix coeff
        Eigen::VectorXd &,
        Eigen::VectorXd& // row order indexing
    ); // pivoting
    static void LUdecompositionTriDiag(
        Eigen::VectorXd &,
        Eigen::VectorXd &,
        Eigen::VectorXd &
    );
};

#endif
