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

#include "basis.h"

using std::clog;
using std::endl;

map<string, string> Basis::basisFile; // baseId -> file
map<string, string> Basis::basisType; // baseId -> type

std::ostream &operator<<(std::ostream &out, Basis &basis)
{
  out << basis.toString();
  return out;
}

//
// Fetch list of available basis sets.
//
// fname - file with basis sets declaration (IN)
//

void Basis::checkBasis(const string &fname)
{
  // Align log to right with dots.
  clog.fill('.');
  clog.width(35);
  clog << std::left << "Loading basis list";

  // Open basis index file.
  std::ifstream f(fname.c_str());

  if (!f.good())
  {
    throw string("error: can not open basis index");
  }

  while (!f.eof())
  {
    string symbol; // basis id
    string file;   // file name, where basis params are stored
    string type;   // basis type (GTO, STO, etc.)

    // Read one basis declaration.
    f >> symbol >> file >> type;

    // Save basis.
    basisFile[symbol] = file;
    basisType[symbol] = type;
  }

  clog << "O.K!" << endl;

  f.close();
};
