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

#include <iomanip>
#include "geometry.h"
#include "atom.h"

Atom::Atom(const Vec3 &r, int Z)
{
  this->r = r;
  this->Z = Z;
};

Atom::Atom(double x, double y, double z, int Z)
{
  this->r = Vec3(x, y, z);
  this->Z = Z;
};

std::ostream &operator<<(std::ostream &out, Atom &atom)
{
  out << std::fixed << std::setw(4) << Geometry::ZToSym(atom.Z) << "   " << atom.r;

  return out;
}

int Atom::getZ()
{
  return Z;
}

Vec3 Atom::getR()
{
  return r;
}
