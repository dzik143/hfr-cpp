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

#ifndef BASIS_H
#define BASIS_H

#include <algebra/algebra.h>
#include <geometry/geometry.h>
#include <psi/psi.h>

class Basis
{
  friend std::ostream &operator<<(std::ostream &, Basis &);

private:

  // Array of atomic orbitals created for
  // *CONCRETE* geometry (molecule).
  vector<Psi> AO;

public:

  virtual void initAO(Geometry &) = 0;
  virtual void computeOverlap(Matrix &) = 0;
  virtual void computeOneElectron(Matrix &, Geometry &) = 0;
  virtual void computeTwoElectron(Tensor4 &) = 0;
  virtual void loadBasis(string &fname) = 0;

  virtual string toString() = 0;
  virtual string showAO() = 0;

  virtual int getDimension() = 0;

  static void checkBasis(const string &fname = "basis.lib");

  static map<string, string> basisFile;
  static map<string, string> basisType;
};

std::ostream &operator<<(std::ostream &out, Basis &baza);

#endif
