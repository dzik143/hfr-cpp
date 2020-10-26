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

#include "vector.h"
#include "matrix.h"

Vector::Vector(int len) : Matrix(len, 1)
{
};

Vector::~Vector()
{
}

double &Vector::operator()(int i)
{
  return data[i];
}

double Vector::operator()(int i) const
{
  return data[i];
};

double &Vector::operator[](int i)
{
  return data[i];
}

double Vector::operator[](int i) const
{
  return data[i];
};

int Vector::len()
{
  // Vector can be horizontal or vertical.
  return std::max(rowsCnt, colsCnt);
}
