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

#include <limits>
#include "hfr.h"

using std::clog;
using std::endl;

HartreeFockRoothan::HartreeFockRoothan()
{
};

void HartreeFockRoothan::solve(Basis &baza, Geometry &geom)
{
  if (geom.getNumberOfElectrons() % 2)
  {
    throw string("error: odd number of electrons");
  }

  ileOrb   = baza.getDimension();             // total number of atomic orbitals
  occupied = geom.getNumberOfElectrons() / 2; // number of occupied orbitals

  C_ = new Matrix(ileOrb, ileOrb);                  // LCAO coefficients
  h_ = new Matrix(ileOrb, ileOrb);                  // one-electron integrals
  G_ = new Tensor4(ileOrb, ileOrb, ileOrb, ileOrb); // two-electrons integrals
  S_ = new Matrix(ileOrb, ileOrb);                  // overlap integrals

  P_ = new Matrix(ileOrb, ileOrb); // bond-order matrix
  F_ = new Matrix(ileOrb, ileOrb); // Fock matrix
  E_ = new Vector(ileOrb);         // orbital energies

  Matrix &C = *C_;
  Matrix &h = *h_;
  Tensor4 &G = *G_;
  Matrix &S = *S_;

  Matrix &P = *P_;
  Matrix &F = *F_;
  Vector &E = *E_;

  // Align log line with dots up to 35 chars.
  clog.fill('.');
  clog.width(35);

  clog << "Calculating molecular integrals";
  baza.computeOverlap(S);
  baza.computeOneElectron(h, geom);
  baza.computeTwoElectron(G);
  clog << "O.K!" << endl;

  Matrix FxS(ileOrb, ileOrb);
  Matrix SFS(ileOrb, ileOrb);
  Matrix Cij_(ileOrb, ileOrb);
  Matrix Cij(ileOrb, ileOrb);

  double Eold, E0, deltaE;
  double THRE;
  double Etotal;

  double NuclearRepulsion;

  int iter, i, j;

  // Init states order array.
  // We start from [0,1,2,3,4,5...] order.
  state.clear();

  for (i = 0; i < ileOrb; i++)
  {
    state.push_back(i);
  }

  // ==========================
  // = Start from zero matrix =
  // ==========================

  for (i = 0; i < ileOrb; i++)
  {
    for (j = 0; j < ileOrb; j++)
    {
      C(i, j) = 0.0;
    }
  }

  sort(state.begin(), state.end(), cmpByE);

  THRE = 1.0e-7; // convergence threshold
  S.power(-0.5); // S = S^(-1/2)

  deltaE = std::numeric_limits<double>::max();
  iter   = 1;

  clog << "Falling into SCF loop..." << endl << endl;

  while (deltaE > THRE)
  {
    setUpP(P, C);          // P = bond-order matrix
    setUpFock(F, h, G, P); // F = Fock-matrix

    // =======================
    // = Solve HFR equation: =
    // = F x C = S x C x E   =
    // =======================

    FxS = F * S;              // FxS = F x S^(-1/2)
    SFS = S * F * S;          // SFS = S^(-1/2) x F x S^(-1/2)
    SFS.diagonalize(E, Cij_); // Diagonalize SFS
    C = S * Cij_;             // Back to original base

    // ==========================
    // = Calculate total energy =
    // ==========================

    // Sort states from the lowest energy.
    sort(state.begin(), state.end(), cmpByE);

    setUpP(P, C);
    setUpFock(F, h, G, P);

    E0 = 0.0;

    for (i = 0; i < ileOrb; i++)
    {
      for (j = 0; j < ileOrb; j++)
      {
        E0 = E0 + P(j, i) * (h(i, j) + F(i, j));
      }
    }

    // Restricted HFR - two electrons per one state.
    E0 = 0.5 * E0;

    // ===============================
    // = Check convergent conditions =
    // ===============================

    Etotal = E0 + geom.nuclearRepulsion();
    deltaE = Eold - Etotal;

    clog.fill('0');
    clog.precision(8);

    clog << "iter=" << iter
         << "   Etotal=" << std::fixed << std::setw(14) << std::left << Etotal;

    if (iter == 1)
    {
      clog << endl;
    }
    else
    {
      clog << " delta=" << std::scientific << deltaE << endl;
    }

    Eold = Etotal;
    iter = iter + 1;
  }

  clog << endl << "SCF convergented." << endl
       << endl
       << "Molecular orbitals:" << endl << C << endl
       << "Orbital energies:" << endl << E;

  delete C_;
  delete h_;
  delete G_;
  delete S_;
  delete P_;
  delete F_;
  delete E_;
}

//
// Build bond-order matrix:
// Prs = Crj * Csj
//

void HartreeFockRoothan::setUpP(Matrix &P, Matrix &C)
{
  int r, s, i, j;
  double suma;

  for (r = 0; r < ileOrb; r++)
  {
    for (s = 0; s < ileOrb; s++)
    {
      suma = 0.0;

      for (i = 0; i < occupied; i++)
      {
        j = state[i];
        suma = suma + C(r, j) * C(s, j);
      }
      P(s, r) = 2.0 * suma;
    }
  }
}

void HartreeFockRoothan::setUpFock(Matrix &F, Matrix &h, Tensor4 &G, Matrix &P)
{
  double suma;

  int r, s, p, q;

  for (r = 0; r < ileOrb; r++)
  {
    for (s = r; s < ileOrb; s++)
    {
      suma = h(r, s);

      for (p = 0; p < ileOrb; p++)
      {
        for (q = 0; q < ileOrb; q++)
        {
          suma = suma + P(q, p) * (G(r, p, s, q) - 0.5 * G(r, p, q, s));
        }
      }

      F(r, s) = suma;
      F(s, r) = suma;
    }
  }
}

Vector *HartreeFockRoothan::E_;

bool HartreeFockRoothan::cmpByE(int i, int j)
{
  return (*E_)(i) < (*E_)(j);
}
