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

#ifndef JOB_H
#define JOB_H

#include <iostream>
#include <sstream>
#include <string>

#include <basis/basis.h>
#include <basis/gto.h>
#include <geometry/geometry.h>
#include <methods/hfr.h>
#include <methods/ssolver.h>

class Job
{

private:

  Geometry geom; // geometria czasteczki
  Basis *baza;   // baza funkcyjna
  Psi *psi;      // funkcja falowa

  std::ifstream *input;

public:

  Job();
  ~Job();

  Job(const string &fname);

  void run();
};

#endif
