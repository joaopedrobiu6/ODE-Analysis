#include <Eigen/Dense>
#include <iostream>
#include "EqSolver.h"
#include "FCmatrixAlgo.h"
#include <math.h>
#include <cmath>

using namespace std;

EqSolver::EqSolver(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &M_, const Eigen::VectorXd &b_)
{
  M = M_;
  b = b_;
  sol1 = b_;
};

const Eigen::VectorXd &EqSolver::GaussSolver(bool pivot)
{
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> M1 = M;
  Eigen::VectorXd b1 = b;

  if(M1.rows()!=b1.rows())
    cout << "\nA Matriz e o Vetor não tem dimensões que permitam resolver o sistema (solução definida)...\n" << endl;

  Eigen::VectorXd &sol = sol1;

  if (pivot == true)
  {
    Eigen::VectorXd piv(M1.rows());

    FCmatrixAlgo::GaussEliminationPivot(M1, b1, piv);
  }
  else if (pivot == false)
  {
    FCmatrixAlgo::GaussElimination(M1, b1);
  }

  // M*sol=b

  sol((M1.rows() - 1)) = (b1(M1.rows() - 1)) / M1((M1.rows() - 1), (M1.cols() - 1));

  for (int i = M1.rows() - 2; i >= 0; i--)
  {
    for (int k = i + 1; k < M1.cols(); k++)
    {
      b1(i) = b1(i) - M1(i, k) * sol(k);
    }
    sol(i) = b1(i) / M1(i, i);
  }

  return sol;
};

const Eigen::VectorXd &EqSolver::LUSolver(bool pivot)
{
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> M1 = M;
  Eigen::VectorXd b1 = b;

  Eigen::VectorXd &sol = sol1;
  Eigen::VectorXd y(b1.rows());
  y.setZero();

  if (pivot == true)
  {
    Eigen::VectorXd piv(M1.rows());
    FCmatrixAlgo::LUdecompositionPivot(M1, b1, piv);
  }
  else if (pivot == false)
  {
    FCmatrixAlgo::LUdecomposition(M1, b1);
  }

  // forward solution (Ly=b)

  y(0) = b1(0);

  for (int i = 1; i < M1.rows(); i++)
  {
    for (int k = i - 1; k >= 0; k--)
    {
      y(i) = y(i) - y(k) * M1(i, k);
    }
    y(i) = b1(i) + y(i);
  }

  cout << "y_1:\n"
       << y << endl;

  // backward solution (Ux=y)

  // loop on rows

  sol((M1.rows() - 1)) = (y(M1.rows() - 1)) / M1((M1.rows() - 1), (M1.cols() - 1));

  for (int i = M1.rows() - 2; i >= 0; i--)
  {
    for (int k = i + 1; k < M1.cols(); k++)
    {
      y(i) = y(i) - M1(i, k) * sol(k);
    }
    sol(i) = y(i) / M1(i, i);
  }

  return sol;
};

void EqSolver::IterativeJacobiSolver(Eigen::VectorXd &IterIn, const int &itmax, double tol)
{
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> M1 = M;
  Eigen::VectorXd b1 = b;

  Eigen::VectorXd aux_sol;

  bool btol = false;
  int it = 0;

  while (!btol && (it++ < itmax))
  {
    aux_sol = IterIn;
    btol = true;
    for (int i = 0; i < M1.rows(); i++)
    {
      IterIn(i) = 0.;
      for (int j = 0; j < M1.rows(); j++)
      {
        if (i != j)
        {
          IterIn(i) = IterIn(i) - M1(i, j) * aux_sol(j);
        }
      }
      IterIn(i) = IterIn(i) + b1(i);
      IterIn(i) = IterIn(i) / M1(i, i);
      // guarantee that all vector entries are converging equally

      if (btol)
      {
        if (fabs(IterIn(i) - aux_sol(i)) < tol)
        {
          btol = true;
        }
        else
        {
          btol = false;
        }
      }
    }
  }
  cout << "\nNúmero de Iterações Realizadas:" << it << endl;
};

void EqSolver::IterativeJacobiSolverRelax(Eigen::VectorXd &IterIn, const int &itmax, double tol, int first, int second)
{
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> M1 = M;
  Eigen::VectorXd b1 = b;

  Eigen::VectorXd aux_sol(b1.rows());

  Eigen::VectorXd x_first(b1.rows());
  Eigen::VectorXd delta_first(b1.rows());
  Eigen::VectorXd delta_second(b1.rows());

  Eigen::VectorXd relax(b1.rows());

  bool btol = false;
  int it = 0;

  while (!btol && (it++ < itmax))
  {
    btol = true;
    aux_sol = IterIn;
    for (int i = 0; i < M1.rows(); i++)
    {
      IterIn(i) = 0;
      for (int j = 0; j < M1.rows(); j++)
      {
        if (i != j)
        {
          IterIn(i) -= M1(i, j) * aux_sol(j);
        }
        IterIn(i) += b1(i);
      }

      if (it <= first + second)
      {
        IterIn(i) /= M1(i, i);
        if (it == first)
        {
          x_first(i) = IterIn(i);
          delta_first(i) = fabs(IterIn(i) - aux_sol(i));
        }
        if (it == first + second)
        {
          delta_second(i) = fabs(IterIn(i) - x_first(i));
          relax(i) = 2 / (1 + sqrt(1 - pow(delta_first(i) / delta_second(i), 1. / (float)second)));
        }
      }
      else
      {
        IterIn(i) *= relax(i) / M1(i, i);
        IterIn(i) += (1 - relax(i)) * aux_sol(i);
      }

      if (btol)
      {
        if (fabs(IterIn(i) - aux_sol(i)) < tol)
        {
          btol = true;
        }
        else
        {
          btol = false;
        }
      }
    }
  }
  cout << "\nNúmero de Iterações Realizadas:" << it << endl;
};

void EqSolver::IterativeGaussSeidelSolver(Eigen::VectorXd &x, int &itmax, double tol)
{
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> M1 = M;
  Eigen::VectorXd b1 = b;

  Eigen::VectorXd aux_x;
  x.setZero();
  aux_x.setZero();
  // linear system of m unknowns Vec x(m); // zero's Vec x_aux(m);     // zero's

  bool btol = false;
  int it = 0.;
  while (!btol && (it++ < itmax))
  {
    aux_x = x;
    for (int i = 0; i < M1.rows(); i++)
    {
      x(i) = 0.;
      for (int j = 0; j < M1.rows(); j++)
        if (i != j)
        {
          x(i) += -M1(i, j) * x(j);
        }

      x(i) += b1(i);
      x(i) /= M1(i, i);

      if (fabs(x(i) - aux_x(i)) < tol)
      {
        btol = true;
      }
      else
      {
        btol = false;
      }
    }
  }
  cout << "\nNúmero de Iterações Realizadas:" << it << endl;
};

