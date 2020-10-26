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

#ifndef TENSOR4_H
#define TENSOR4_H

class Tensor4
{

protected:

  double *data;

  int iMax_;
  int jMax_;
  int kMax_;
  int lMax_;

public:

  Tensor4(int = 0, int = 0, int = 0, int = 0);

  ~Tensor4();

  // Fortran-like Access operators A(i,j,k,l)
  double &operator()(int, int, int, int);
  double operator()(int, int, int, int) const;

  Tensor4 &operator=(const Tensor4 &);
  Tensor4 operator+(const Tensor4 &);
  Tensor4 operator-(const Tensor4 &);
  Tensor4 operator*(const Tensor4 &);

  int iMax();
  int jMax();
  int kMax();
  int lMax();

  bool empty();
};

#endif
