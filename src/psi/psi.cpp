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

#include "psi.h"

map<string, int> Psi::typToL; // s,p,d... -> 0,1,2...

Psi::Psi()
{
  typToL["S"] = 0;
  typToL["P"] = 1;
  typToL["D"] = 2;
  typToL["F"] = 3;
  typToL["G"] = 4;
} //Psi

std::ostream &operator<<(std::ostream &out, Psi &psi)
{
  out << psi.toString();
  return out;
}

//
// Contracted gauss.
// Psi = W*[c1*Exp(-a1*r1) + c2*Exp(-a2*r2) + ...]
//

GTO::GTO(vector<double> &a,vector<double> &c,
         const string &typ, const string &orient, const Vec3 &r)
{
  this->a = a;                     // coefficient in exponent (alpha)
  this->c = c;                     // coefficient before exponent
  this->lTotal = Psi::typToL[typ]; // s,p,d... -> 0,1,2..
  this->r = r;                     // center coordinates
}

GTO::GTO(const string &typ, const string &orient,
         const Vec3 &r)
{
  this->lTotal = Psi::typToL[typ];
  this->r = r;
}

GTO::GTO(vector<double> &a, vector<double> &c,
         int lx, int ly, int lz, const Vec3 &r)
{
  this->a = a;                 // exponent coeffcient (alpha)
  this->c = c;                 // center coordinates
  this->lTotal = lx + ly + lz; // s,p,d,... -> 0,1,2...
  this->r = r;                 // center coordinates
}

//
// Add one primitive GTO to contracted set: c * exp(ax)
//

void GTO::addPrimitive(double c, double a)
{
  this->c.push_back(c);
  this->a.push_back(a);
}

//
// Dump contracted GTO list to the output stream as
// human-readable text.
//

string GTO::toString()
{
  std::ostringstream ss;

  ss << "x^" << lx << " y^" << ly << " z^" << lz << " x ";

  for (unsigned int i = 0; i < c.size(); i++)
  {
    ss << "(" << c[i] << "," << a[i] << ") ";
  }

  ss << "}";

  return ss.str();
}

//
// Move GTO center by vector.
//

void GTO::przesun(Vec3 &dR)
{
  r = r + dR;
}

//
// Getters.
//

int GTO::getL() { return lTotal; }

Vec3 GTO::getR() { return r; }

vector<double> &GTO::getC() { return c; }
vector<double> &GTO::getA() { return a; }
