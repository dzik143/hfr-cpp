/******************************************************************************/
/*                                                                            */
/* Copyright (C) 2009 Sylwester Wysocki <sw143@wp.pl>                         */
/*                                                                            */
/* Source code available at: https://github.com/dzik143/hfr-cpp               */
/*                                                                            */
/* This program is free software: you can redistribute it and/or modify       */
/* it under the terms of the GNU General Public License as published by       */
/* the Free Software Foundation, either version 3 of the License, or          */
/* (at your option) any later version.                                        */
/*                                                                            */
/* This program is distributed in the hope that it will be useful,            */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of             */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              */
/* GNU General Public License for more details.                               */
/*                                                                            */
/* You should have received a copy of the GNU General Public License          */
/* along with this program. If not, see <http://www.gnu.org/licenses/>        */
/*                                                                            */
/******************************************************************************/

#include <cmath>
#include <iomanip>

#include "matrix.h"
#include "error.h"
#include "vector.h"

//
// Init code: constructors
//

Matrix::Matrix(int rowsCnt, int colsCnt)
{
  data = NULL;

  if ((rowsCnt != 0) && (colsCnt != 0))
  {
    data = new double[rowsCnt * colsCnt];

    if (data == NULL)
    {
      throw Error("Out of memory", "ALGEBRA", "Matrix::Matrix");
    }
  };

  this->rowsCnt = rowsCnt;
  this->colsCnt = colsCnt;
}

Matrix::Matrix(int rowsCnt, int colsCnt, const double **initialData)
{
  data = NULL;

  if ((rowsCnt != 0) && (colsCnt != 0))
  {
    data = new double[rowsCnt * colsCnt];

    if (data == NULL)
    {
      throw Error("Out of memory", "ALGEBRA", "Matrix::Matrix");
    }
  }

  for (int i = 0; i < rowsCnt; i++)
  {
    for (int j = 0; j < rowsCnt; j++)
    {
      (*this)(i, j) = initialData[i][j];
    }
  }

  this->rowsCnt = rowsCnt;
  this->colsCnt = colsCnt;
}

Matrix::Matrix(const Matrix &B)
{
  if (this == &B)
  {
    return;
  }

  if (colsCnt * rowsCnt < B.colsCnt * B.rowsCnt)
  {
    if (data != NULL)
    {
      delete data;
    }
    data = new double[B.colsCnt * B.rowsCnt];

    if (data == NULL)
    {
      throw Error("Out of memory", "ALGEBRA", "Matrix::operator=");
    }
  }

  (*this).rowsCnt = B.rowsCnt;
  (*this).colsCnt = B.colsCnt;

  for (int i = 0; i < rowsCnt; i++)
  {
    for (int j = 0; j < colsCnt; j++)
    {
      (*this)(i, j) = B(i, j);
    }
  }
}

//
// Clean up: destructor
//

Matrix::~Matrix()
{
  if (data != NULL)
  {
    delete data;
  }
}

//
// Assign operator
// A = B
//

Matrix &Matrix::operator=(const Matrix &B)
{
  if (this == &B)
  {
    return *this;
  }

  if (colsCnt * rowsCnt < B.colsCnt * B.rowsCnt)
  {
    if (data != NULL)
    {
      delete data;
    }

    data = new double[B.colsCnt * B.rowsCnt];

    if (data == NULL)
    {
      throw Error("Out of memory", "ALGEBRA", "Matrix::operator=");
    }
  }

  (*this).rowsCnt = B.rowsCnt;
  (*this).colsCnt = B.colsCnt;

  for (int i = 0; i < rowsCnt; i++)
  {
    for (int j = 0; j < colsCnt; j++)
    {
      (*this)(i, j) = B(i, j);
    }
  }

  return *this;
}

//
// Getters.
//

int  Matrix::szer()  { return colsCnt; }
int  Matrix::wys()   { return rowsCnt; }
bool Matrix::empty() { return (data == NULL); }

//
// Fortran like access operator:
//
// A(i,j) = ...
// ...    = A(i,j)
//
//

double &Matrix::operator()(int w, int k)
{
  if ((w >= rowsCnt) || (k >= colsCnt))
  {
    throw Error("Out of range", "ALGEBRA", "Matrix::operator()");
  }

  return data[colsCnt * w + k];
}

double Matrix::operator()(int w, int k) const
{
  if ((w >= rowsCnt) || (k >= colsCnt))
  {
    throw Error("Out of range", "ALGEBRA", "Matrix::operator()");
  }
  return data[colsCnt * w + k];
}

// -----------------------------------------------------------------------------
//                                 I/O related
// -----------------------------------------------------------------------------

//
// Read matrix from the input stream e.g. from the file.
//

std::istream &operator>>(std::istream &inp, Matrix &arg)
{
  if (arg.empty())
  {
    throw Error("Empty target matrix", "ALGEBRA", "Matrix::operator>>");
  }

  for (int j = 0; j < arg.wys(); j++)
  {
    for (int i = 0; i < arg.szer(); i++)
    {
      inp >> arg(j, i);
    }
  }

  return inp;
};

//
// Write matrix to the output stream e.g. to the file.
//

std::ostream &operator<<(std::ostream &out, Matrix &arg)
{
  out.precision(3);
  out.fill(' ');

  for (int j = 0; j < arg.wys(); j++)
  {
    for (int i = 0; i < arg.szer(); i++)
    {
      out << fixed << std::setw(8) << std::right << arg(j, i);
    }

    out << "\n";
  }
  return out;
}

// -----------------------------------------------------------------------------
//                           Basic algebra operations
// -----------------------------------------------------------------------------

//
// Matrix add: C = A + B
//

Matrix Matrix::operator+(const Matrix &B)
{
  Matrix C(rowsCnt, colsCnt);

  if ((colsCnt != B.colsCnt) || (rowsCnt != B.rowsCnt))
  {
    throw Error("Bad dimensions", "ALGEBRA", "Matrix::operator+");
  }

  for (int i = 0; i < rowsCnt; i++)
  {
    for (int j = 0; j < colsCnt; j++)
    {
      C(i, j) = (*this)(i, j) + B(i, j);
    }
  }

  return C;
}

//
// Matrix subtract: C = A - B
//

Matrix Matrix::operator-(const Matrix &B)
{
  Matrix C(rowsCnt, colsCnt);

  if ((colsCnt != B.colsCnt) || (rowsCnt != B.rowsCnt))
  {
    throw Error("Bad dimensions", "ALGEBRA", "Matrix::operator-");
  }

  for (int i = 0; i < rowsCnt; i++)
  {
    for (int j = 0; j < colsCnt; j++)
    {
      C(i, j) = (*this)(i, j) - B(i, j);
    }
  }

  return C;
}

//
// General matrix multiplication: C = A x B
//

Matrix Matrix::operator*(const Matrix &B)
{
  Matrix C(rowsCnt, B.colsCnt);

  if (colsCnt != B.rowsCnt)
  {
    throw Error("Bad dimensions", "ALGEBRA", "Matrix::operator*");
  }

  for (int i = 0; i < rowsCnt; i++)
  {
    for (int j = 0; j < B.colsCnt; j++)
    {
      double suma = 0.0;

      for (int k = 0; k < colsCnt; k++)
      {
        // Cij = Aik * Bkj
        suma = suma + (*this)(i, k) * B(k, j);
      }

      C(i, j) = suma;
    }
  }

  C.rowsCnt = (*this).rowsCnt;
  C.colsCnt = B.colsCnt;

  return C;
}

//
// Matrix by scalar mutliplication: C = A * s
//

Matrix Matrix::operator*(const double s)
{
  Matrix C(rowsCnt, colsCnt);

  for (int i = 0; i < rowsCnt; i++)
  {
    for (int j = 0; j < colsCnt; j++)
    {
      C(i, j) = (*this)(i, j) * s;
    }
  }

  return C;
}

// -----------------------------------------------------------------------------
//                                   Other
// -----------------------------------------------------------------------------

//
// Helper function to get argumnent sign.
//

double Matrix::sgn(double x)
{
  if (x < 0.0)
  {
    return -1.0;
  }

  return 1.0;
}

//
// Diagonalize matrix using Jacobi method:
// https://en.wikipedia.org/wiki/Jacobi_method
//
// WARNING #1: Method destroys values stored in the matrix.
// WARNING #2: Result eigenvalues are *NOT* sorted.
//
// eVal - eigenvalues (OUT)
// eVec - eigenvectors (OUT)
//

void Matrix::diagonalize(Matrix &eVal, Matrix &eVec)
{
  int dim = rowsCnt;
  int p, q;
  int i, j;

  Matrix &T = (*this);
  Matrix S(dim, dim);

  // Helper variables.
  double u, v, lambda;
  double sinx, cosx, sinxcos, sin2, cos2;
  double Tpq, Tqq, Tpp, Tpj, Tip, Tmax, Sip, Siq;

  double thre;
  double teta, tau, te;

  int iter = 0;

  const int maxIter = 5000;

  // ============
  // = S0 = [1] =
  // ============

  // Load identity matrix at the beggining.
  for (i = 0; i < dim; i++)
  {
    for (j = 0; j < dim; j++)
    {
      if (i == j)
      {
        S(i, j) = 1.0;
      }
      else
      {
        S(i, j) = 0.0;
      }
    }
  }

  Tmax = 1000.0;
  thre = 1e-5;

  //
  // Go on as long as maximum is out-of the diagonal.
  //

  while ((Tmax > thre) && (iter < maxIter))
  {
    // ====================
    // = szukaj (p,q) max =
    // ====================

    Tmax = 0.0;

    for (i = 0; i < dim; i++)
    {
      for (j = i + 1; j < dim; j++)
      {
        // Skip the diagonal.
        if (abs(T(i, j)) > Tmax)
        {
          // New maximum value.
          Tmax = abs(T(i, j));
          p = i;
          q = j;
        }
      }
    }

    if (Tmax == 0.0)
    {
      break;
    }

    // =================
    // = Helper values =
    // =================

    teta = (T(q, q) - T(p, p)) / (2.0 * T(p, q));
    te   = sgn(teta) / (abs(teta) + sqrt(teta * teta + 1.0));
    cosx = 1.0 / sqrt(te * te + 1.0);
    sinx = te * cosx;
    tau  = sinx / (1. + cosx);

    sinxcos = sinx * cosx;
    sin2    = sinx * sinx;
    cos2    = cosx * cosx;

    Tpp = cos2 * T(p, p) + sin2 * T(q, q) - 2. * sinxcos * T(p, q);
    Tqq = sin2 * T(p, p) + cos2 * T(q, q) + 2. * sinxcos * T(p, q);
    Tpq = (cos2 - sin2) * T(p, q) + sinxcos * (T(p, p) - T(q, q));

    // ==============================
    // = Go through row (P,j),(Q,j) =
    // ==============================

    for (j = 0; j < dim; j++)
    {
      Tpj = T(p, j);
      T(p, j) = Tpj * cosx - T(q, j) * sinx;
      T(q, j) = Tpj * sinx + T(q, j) * cosx;
    }

    // =================================
    // = Go through column (i,P),(i,Q) =
    // =================================

    for (i = 0; i < dim; i++)
    {
      Tip = T(i, p);
      T(i, p) = cosx * Tip - sinx * T(i, q);
      T(i, q) = sinx * Tip + cosx * T(i, q);

      Sip = S(i, p);
      Siq = S(i, q);
      S(i, p) = Sip * cosx - Siq * sinx; // eigen
      S(i, q) = Sip * sinx + Siq * cosx; // vectors
    }

    T(p, q) = Tpq;
    T(q, p) = Tpq;
    T(p, p) = Tpp;
    T(q, q) = Tqq;

    iter++;
  }

  if (iter >= maxIter)
  {
    throw Error("Too many iters", "ALGEBRA", "Matrix::diagonalize");
  }

  // ====================================
  // = Set up result eigenvectors array =
  // ====================================

  for (i = 0; i < dim; i++)
  {
    eVal[i] = T(i, i);
  }

  eVec = S;
}

//
// Raise matrix to given power: C = A^n
//

void Matrix::power(double n)
{
  int dim = rowsCnt;
  int i, j;

  Matrix &T = (*this);
  Vector eValue(dim);
  Matrix eVec(dim, dim);

  T.diagonalize(eValue, eVec);

  for (i = 0; i < dim; i++)
  {
    if (T(i, i) <= 0.0)
    {
      throw Error("Singular matrix", "ALGEBRA", "Matrix::power");
    }

    T(i, i) = pow(T(i, i), n);
  }

  T = eVec * T * eVec.transpone();
}

//
// Transpone matrix: A(i,j) = A(j,i)
//

Matrix Matrix::transpone()
{
  int dim = rowsCnt;

  Matrix wynik(dim, dim);

  for (int i = 0; i < dim; i++)
  {
    for (int j = 0; j < dim; j++)
    {
      wynik(i, j) = (*this)(j, i);
    }
  }

  return wynik;
}

//
// Access to raw underlying buffer.
//

double &Matrix::operator[](int i)
{
  return data[i];
}

double Matrix::operator[](int i) const
{
  return data[i];
};

//
// Invert matrix using Gauss-Jordan method: B = A^-1
// https://www.mathsisfun.com/algebra/matrix-inverse-row-operations-gauss-jordan.html
//

Matrix Matrix::invert()
{
  int dim = rowsCnt;

  Matrix Ai(dim, dim);
  Matrix A(dim, dim);

  A = (*this);

  // ==============================
  // = Start from identity matrix =
  // ==============================

  for (int i = 0; i < dim; i++)
  {
    for (int j = 0; j < dim; j++)
    {
      if (i == j)
      {
        Ai(i, j) = 1.0;
      }
      else
      {
        Ai(i, j) = 0.0;
      }
    }
  }

  for (int j = 0; j < dim; j++)
  {
    // =======================================
    // = If zero on diagonal, then swap rows =
    // =======================================

    if (A(j, j) == 0.0)
    {
      if (j == dim)
      {
        throw Error("Singular matrix", "ALGEBRA", "Matrix::invert");
      }

      int k = j;

      // Search a row *WITHOUT* zero value.
      while ((A(k, j) == 0.0) && (k < dim))
      {
        k++;
      }

      // ======================
      // = swap rows: j <=> l =
      // ======================

      for (int l = 0; l < dim; l++)
      {
        std::swap(A(j, l), A(k, l));
        std::swap(Ai(j, l), Ai(k, l));
      }
    }

    // ====================================
    // = If still zero, then swap columns =
    // ====================================

    if (A(j, j) == 0.0)
    {
      int k = j;

      while ((A(j, k) == 0.0) && (k < dim))
      {
        k++;
      }

      // =========================
      // = swap columns: j <=> l =
      // =========================

      for (int l = 0; l < dim; l++)
      {
        std::swap(A(l, j), A(l, k));
        std::swap(Ai(l, j), Ai(l, k));
      }
    }

    // ==========================================
    // = Scale a row to get 1.0 on the diagonal =
    // ==========================================

    double scale = 1.0 / A(j, j);

    for (int l = 0; l < dim; l++)
    {
      A(j, l)  = A(j, l)  * scale;
      Ai(j, l) = Ai(j, l) * scale;
    }

    for (int i = 0; i < dim; i++)
    {
      if (i != j)
      {
        double scale = A(i, j) / A(j, j);

        for (int l = 0; l < dim; l++)
        {
          A(i, l)  = A(i, l)  - A(j, l)  * scale;
          Ai(i, l) = Ai(i, l) - Ai(j, l) * scale;
        }
      }
    }
  }

  return Ai;
}

//
// Load identity matrix:
// Aii = 1.0
// Aij = 0.0
//

void Matrix::loadIdentity()
{
  int i, j;

  for (i = 0; i < rowsCnt; i++)
  {
    for (j = 0; j < colsCnt; j++)
    {
      if (i == j)
      {
        (*this)(i, j) = 1.0;
      }
      else
      {
        (*this)(i, j) = 0.0;
      }
    }
  }
}
