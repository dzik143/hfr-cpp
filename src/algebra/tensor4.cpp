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

#include "tensor4.h"
#include "error.h"

//
// Init: constructor
//

Tensor4::Tensor4(int iMax, int jMax, int kMax, int lMax)
{
  // Allocate buffer.
  data = new double[iMax * jMax * kMax * lMax];

  if (data == NULL)
  {
    throw Error("Out of memory", "ALGEBRA", "Tensor4::Tensor4");
  }

  // Save dimensions.
  this->iMax_ = iMax;
  this->jMax_ = kMax;
  this->kMax_ = lMax;
  this->lMax_ = lMax;
}

//
// Clean up: destructor
//

Tensor4::~Tensor4()
{
  if (data != NULL)
  {
    delete[] data;
  }
}

//
// Fortran-like access operator:
// A(i,j,k,l) = ...
// ...        = A(i,j,k,l)
//

double &Tensor4::operator()(int i, int j, int k, int l)
{
  return data[i + (iMax_ * (j + jMax_ * (k + l * kMax_)))];
}

double Tensor4::operator()(int i, int j, int k, int l) const
{
  return data[i + (iMax_ * (j + jMax_ * (k + l * kMax_)))];
}

//
// Getters.
//

int  Tensor4::iMax()  { return iMax_; }
int  Tensor4::jMax()  { return jMax_; }
int  Tensor4::kMax()  { return kMax_; }
int  Tensor4::lMax()  { return lMax_; }
bool Tensor4::empty() { return (data == NULL); }
