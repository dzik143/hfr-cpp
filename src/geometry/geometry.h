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

#ifndef GEOM_H
#define GEOM_H

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "vec3.h"
#include "atom.h"

using std::vector;
using std::map;
using std::string;

class Geometry
{
  friend std::ostream &operator<<(std::ostream &, Geometry &);

private:

  vector<Atom> atomy; // lista atomow

  static map<string, int> symToZMap; // symbol -> liczba. at.
  static map<int, string> ZToSymMap; // symbol <- liczba. at.

  int numberOfElectrons; // ile elektronow w ukladzie

public:

  // Calculate nuclei-nuclei repulsion energy.
  double nuclearRepulsion();

  void loadGeom(std::istream &inp);
  void loadGeom(std::string &fname);

  int ileAtom();                // zwraca l. atomow

  void addAtom(Vec3 &, int);    // dodaje atom do ukladu
  void addAtom(Vec3 &, string); // dodaje atom do ukladu
  void addAtom(double, double, double, int);
  void addAtom(double, double, double, string);
  void eraseAtom(int); // usuwa atom z ukladu

  int getNumberOfElectrons();

  Atom &operator[](int);

  static void loadSymbols(string fname = "symbols.dat");
  static int symToZ(const string sym);
  static string ZToSym(int Z);
};

std::ostream &operator<<(std::ostream &, Geometry &);

#endif
