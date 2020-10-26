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
#include "vec3.h"

Vec3::Vec3(double x, double y, double z)
{
  this->x = x;
  this->y = y;
  this->z = z;
};

std::ostream &operator<<(std::ostream &out, Vec3 &v)
{
  out << std::right << v.x << " " << v.y << " " << v.z;

  return out;
};

//
// Add vectors: C = A + B
//

Vec3 Vec3::operator+(Vec3 &v)
{
  Vec3 wynik;
  wynik.x = this->x + v.x;
  wynik.y = this->y + v.y;
  wynik.z = this->z + v.z;

  return wynik;
}

//
// Calculate distance between two points:
// d = sqrt(dx^2 + dy^2 + dz^2)

double Vec3::dist(Vec3 &v2)
{
  double dx = this->x - v2.x;
  double dy = this->y - v2.y;
  double dz = this->z - v2.z;

  return sqrt(dx * dx + dy * dy + dz * dz);
}

//
// Calculate squared distance between two points:
// d^2 = dx^2 + dy^2 + dz^2
//

double Vec3::dist2(Vec3 &v2)
{
  double dx = this->x - v2.x;
  double dy = this->y - v2.y;
  double dz = this->z - v2.z;

  return dx * dx + dy * dy + dz * dz;
}
