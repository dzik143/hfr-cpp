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

#ifndef ATOM_H
#define ATOM_H

#include <ostream>
#include "vec3.h"

// Forward declaration.
class Geometry;

class Atom
{
  friend class Geometry;
  friend std::ostream &operator<<(std::ostream &, Atom &);

private:

  Vec3 r; // Nuclei center
  int Z;  // Z-number (number of protons in the nuclei)

public:

  Atom(const Vec3 &r, int Z);
  Atom(double x, double y, double z, int Z);

  int getZ();

  Vec3 getR();

  inline double x() { return r.x; }
  inline double y() { return r.y; }
  inline double z() { return r.z; }
};

std::ostream &operator<<(std::ostream &, Atom &);

#endif
