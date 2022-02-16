#include "FCmatrixAlgo.h"
#include <iostream>

using namespace std;

void FCmatrixAlgo::GaussElimination(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &A, Eigen::VectorXd &b)
{
    float lambda;
    for (int i = 0; i < A.cols(); i++)
    {
        for (int j = i + 1; j <= A.rows() - 1; j++)
        {
            lambda = A(j, i) / A(i, i);
            A.row(j) = A.row(j) - lambda * (A.row(i));
            b[j] = b[j] - lambda * (b[i]);

            if (A(j, i) < 1e-05)
            {
                A(j, i) = 0;
            }
        }
    }
};

void FCmatrixAlgo::GaussEliminationPivot(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &A, Eigen::VectorXd &b, Eigen::VectorXd &piv)
{
    for (int a = 0; a < A.rows(); a++)
    {
        if (fabs(A.row(a).maxCoeff()) >= fabs(A.row(a).minCoeff()))
        {
            piv(a, 0) = A.row(a).maxCoeff();
        }
        else
        {
            piv(a, 0) = fabs(A.row(a).minCoeff());
        }
    }

    double scale;
    for (int i = 0; i < A.rows() - 1; i++)
    {

        for (int k = i + 1; k < A.rows(); k++)
        {
            if (fabs(A(i, i) / piv(i, 0)) < fabs(A(k, i) / piv(k, 0)))
            {
                A.row(i).swap(A.row(k));
                b.row(i).swap(b.row(k));
            }
        }
        for (int j = i + 1; j < A.rows(); j++)
        {
            scale = -A(j, i) / A(i, i);
            A.row(j) += scale * A.row(i);
            b(j) += scale * b(i);

            if (A(j, i) < 1e-05)
            {
                A(j, i) = 0;
            }
        }
    }
};

void FCmatrixAlgo::LUdecomposition(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &A, Eigen::VectorXd &b)
{
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> B(A.rows(), A.cols());

    float lambda;
    for (int i = 0; i < A.cols(); i++)
    {
        for (int j = i + 1; j <= A.rows() - 1; j++)
        {
            lambda = A(j, i) / A(i, i);
            A.row(j) = A.row(j) - lambda * (A.row(i));
            B(j, i) = lambda;
        }
    }
    A = A + B;
};

void FCmatrixAlgo::LUdecompositionPivot(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &A, Eigen::VectorXd &b, Eigen::VectorXd &piv)
{
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> B(A.rows(), A.cols());

    for (int a = 0; a < A.rows(); a++)
    {
        if (fabs(A.row(a).maxCoeff()) >= fabs(A.row(a).minCoeff()))
        {
            piv(a, 0) = A.row(a).maxCoeff();
        }
        else
            piv(a, 0) = fabs(A.row(a).minCoeff());
    }

    for (int i = 0; i < A.rows() - 1; i++)
    {
        for (int k = i + 1; k < A.rows(); k++)
        {
            if (fabs(A(i, i) / piv(i)) < fabs(A(k, i) / piv(k)))
            {
                A.row(i).swap(A.row(k));
                b.row(i).swap(b.row(k));
            }
        }
    }

    float lambda;
    for (int i = 0; i < A.cols(); i++)
    {
        for (int j = i + 1; j <= A.rows() - 1; j++)
        {
            lambda = A(j, i) / A(i, i);
            A.row(j) = A.row(j) - lambda * (A.row(i));
            B(j, i) = lambda;
        }
    }
    A = A + B;
};

void FCmatrixAlgo::LUdecompositionTriDiag(Eigen::VectorXd &a, Eigen::VectorXd &b, Eigen::VectorXd &c)
{
    for (int i = 1; i < b.rows(); i++) {
        a(i-1) = a(i-1)/b(i-1);
        b(i) -= a(i-1)*c(i-1);
    }
};
