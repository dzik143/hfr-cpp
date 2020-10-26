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

#ifndef VEC3_H
#define VEC3_H

#include <ostream>

class Vec3
{
  friend std::ostream &operator<<(std::ostream &, Vec3 &);

public:

  double x, y, z;

  Vec3(double x = 0.0, double y = 0.0, double z = 0.0);
  Vec3 operator+(Vec3 &);

  double dist(Vec3 &);
  double dist2(Vec3 &);
};

std::ostream &operator<<(std::ostream &, Vec3 &);

#endif
