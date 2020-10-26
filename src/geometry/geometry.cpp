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

#include "geometry.h"

using std::clog;
using std::endl;

map<string, int> Geometry::symToZMap; // elementSymbol -> Z-number
map<int, string> Geometry::ZToSymMap; // elementSymbol <- Z-number

//
// Load geometry in format:
//   symbol1 x y z
//   symbol2 x y z
//   symbol3 x y z
//   ...
//   END_GEOM
//
// inp - input stream to read geometry in text format (IN)
//

void Geometry::loadGeom(std::istream &inp)
{
  string sym;
  double x, y, z;
  int Z;

  // =================
  // = Load geometry =
  // =================

  clog.fill('.');
  clog.width(35);
  clog << "Loading geometry";

  atomy.clear();
  numberOfElectrons = 0;

  // First atom.
  inp >> std::ws >> sym;

  // Read atoms line-by-line until END_GEOM reached.
  while (sym != "END_GEOM")
  {
    // Read next atom.
    inp >> x >> y >> z;
    Z = Geometry::symToZ(sym);

    if (Z == 0)
    {
      throw string("error: invalid geometry");
    }

    // Track total number of electron in zero-charge system.
    numberOfElectrons += Z;

    // Convert angstrems to a.u.
    x = x / 0.5292;
    y = y / 0.5292;
    z = z / 0.5292;

    // Add atom to the system.
    atomy.push_back(Atom(x, y, z, Z));

    // Fetch next atom.
    inp >> std::ws >> sym;
  }

  // ================
  // = Log geometry =
  // ================

  clog.precision(5);

  clog << "O.K!"
       << endl
       << endl
       << " INPUT GEOMETRY [au]"            << endl
       << "-------------------------------" << endl
       << "  Atom      X      Y      Z"     << endl
       << *this
       << "-------------------------------"
       << endl
       << endl;
};

//
// Load chemical elements symbols from the file in format:
//   <Z-number> <symbol>
//
// fname - path to input file, defaulted to "symbols.dat"
//         if skipped (IN/OPT)
//

void Geometry::loadSymbols(string fname)
{
  std::ifstream finp;

  clog.fill('.');                    // dopelnienie kropkami
  clog.width(35);                    // do 35 znakow
  clog << "Loading element symbols"; //

  finp.open(fname.c_str());

  while (!finp.eof())
  {
    string sym;
    int Z;

    finp >> Z >> sym;             // czytaj l.at. i symbol
    Geometry::symToZMap[sym] = Z; // Zapisz mape z symbolami
    Geometry::ZToSymMap[Z] = sym; //
  };
  finp.close();

  clog << "O.K!" << endl;
};

//
// Convert element symbol to Z-number.
// Example: "He" -> 2 (helium).
//

int Geometry::symToZ(const string sym)
{
  return Geometry::symToZMap[sym];
}

//
// Convert get symbol of element with given Z-number.
// Example: 2 -> He (helium).
//

string Geometry::ZToSym(int Z)
{
  return Geometry::ZToSymMap[Z];
}

//
// Dump geometry to the output stream as human-readable text.
//

std::ostream &operator<<(std::ostream &out, Geometry &geom)
{
  out.fill(' ');
  out << std::right;

  vector<Atom>::iterator it;

  for (it = geom.atomy.begin(); it != geom.atomy.end(); it++)
  {
    out << *it << endl;
  }

  return out;
};

//
// Access to the i-th atom.
//

Atom &Geometry::operator[](int i)
{
  return atomy[i];
}

int Geometry::ileAtom()
{
  return atomy.size();
}

//
// Add one atom to the system.
//
// r - atom center (x,y,z) (IN),
// Z - Z-number of atom e.g. 2 for Helium.
//

void Geometry::addAtom(Vec3 &r, int Z)
{
  atomy.push_back(Atom(r, Z));
  numberOfElectrons += Z;
};

//
// Add new atom to the system.
//
// r   - atom center (x,y,z) (IN),
// sym - element symbol e.g. "He" for helium.
//

void Geometry::addAtom(Vec3 &r, string sym)
{
  atomy.push_back(Atom(r, symToZ(sym)));
  numberOfElectrons += symToZ(sym);
};

//
// Add one atom to the system.
//
// x,y,z - atom center (IN),
// Z     - Z-number of atom e.g. 2 for Helium.
//

void Geometry::addAtom(double x, double y, double z, int Z)
{
  atomy.push_back(Atom(x, y, z, Z));
  numberOfElectrons += Z;
};

//
// Add new atom to the system.
//
// x,y,z - atom center (IN),
// sym   - element symbol e.g. "He" for helium.
//

void Geometry::addAtom(double x, double y, double z, string sym)
{
  atomy.push_back(Atom(x, y, z, symToZ(sym)));
  numberOfElectrons += symToZ(sym);
};

//
// Removes i-th atom from the system.
//

void Geometry::eraseAtom(int i)
{
  if (i >= atomy.size())
  {
    throw string("error: no such atom");
  }

  numberOfElectrons -= atomy[i].getZ();
  atomy.erase(atomy.begin() + i);
};

//
// Calculate sum of nuclear repulsion energy between each nuclei pairs:
//
//       Z1 * Z2
// E12 = -------
//         r12
//

double Geometry::nuclearRepulsion()
{
  double Ejj;
  int i, j;
  int Za, Zb;
  double dx, dy, dz;
  double Rab;

  const double k = 0.52917693824475;

  Ejj = 0.0;

  // jedz po parach jader atomowych
  for (i = 0; i < ileAtom(); i++)
  {
    for (j = i + 1; j < ileAtom(); j++)
    {
      dx = atomy[i].x() - atomy[j].x();
      dy = atomy[i].y() - atomy[j].y();
      dz = atomy[i].z() - atomy[j].z();
      Rab = sqrt(dx * dx + dy * dy + dz * dz);

      Za = atomy[i].Z;
      Zb = atomy[j].Z;

      Ejj = Ejj + Za * Zb / Rab;
    }
  }

  return Ejj;
}

int Geometry::getNumberOfElectrons()
{
  return numberOfElectrons;
}
