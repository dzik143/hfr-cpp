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

#ifndef HFR_H
#define HFR_H

#include <algebra/algebra.h>
#include <basis/basis.h>
#include <geometry/geometry.h>
#include <methods/ssolver.h>
#include <psi/psi.h>

#include <algorithm>
#include <vector>

class HartreeFockRoothan : public SchroedingerSolver
{

public:

  HartreeFockRoothan();

  void solve(Basis &, Geometry &);

private:

  // LCAO coefficients (eigenvectors go here).
  Matrix *C_;

  // Intermediate integral matrices:
  // h - one electron integrals, hpq   = <p|T+V|q>
  // G - two-electron integrals, Grpqs = <pq|G|rs> = <pq| 1/r12 | rs>
  // S - overlap integrals, Spq = <p|q>

  Matrix *S_;  // Overlap
  Matrix *h_;  // one-electron
  Tensor4 *G_; // Two-electrons

  // Intermediate matrixes for FC=SCE equation.
  Matrix *P_; // bond-order matrix
  Matrix *F_; // Fock matrix

  // Total and occupied orbitals.
  int ileOrb;
  int occupied;

  // Orbital energies (eigenvalues).
  static Vector *E_;

  // States order from the lowest energy.
  vector<int> state;

  // Bond-order matrix.
  void setUpP(Matrix &, Matrix &);

  // Fock matrix.
  void setUpFock(Matrix &, Matrix &, Tensor4 &, Matrix &);

  // To sort orbitals (eigenvectors) by energy (eigenvalues).
  static bool cmpByE(int, int);
};

#endif
