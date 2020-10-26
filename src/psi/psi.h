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

#ifndef PSI_H
#define PSI_H

#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <geometry/geometry.h>

//
// Abstract orbital function.
//

class Psi
{
  friend std::ostream &operator<<(std::ostream &, Psi &);

protected:

  static std::map<std::string, int> typToL; // s,p,d... -> 0,1,2...

public:

  Psi();

  virtual string toString() = 0;
  virtual void przesun(Vec3 &) = 0;

  virtual int getL() = 0;
  virtual Vec3 getR() = 0;
};

std::ostream &operator<<(std::ostream &out, Psi &psi);

//
// Primitive GTO function:
// psi = W(x,y,z) * exp(-ar)
//

class PGTO : public Psi
{

private:

  Vec3 r;   // center coordinates
  double a; // alpha coefficient (under exponent)
  string l; // poboczna l.k
  string m; // magnetyczna l.k

public:

  PGTO(const Vec3 &r = Vec3(0, 0, 0), double a = 1.0,
       const string &l = "s", const string &m = "");
};

//
// Contracted GTO function (waged sum of primitive gauss with
// common center):
// Psi = W*[c1*Exp(-a1*r1) + c2*Exp(-a2*r2) + ...]
//

class GTO : public Psi
{

private:

  Vec3 r;           // center coordinates
  int lTotal;       // lx+ly+lz
  int lx, ly, lz;   // powers at xa,ya,za
  vector<double> a; // alpha coefficients (under exponents)
  vector<double> c; // contraction coefficients

public:

  GTO(vector<double> &a, vector<double> &c,
      const string &typ = "s", const string &orient = "",
      const Vec3 &r = Vec3(0, 0, 0));

  GTO(const string &typ, const string &orient = "",
      const Vec3 &r = Vec3(0, 0, 0));

  GTO(vector<double> &a, vector<double> &c,
      int lx, int ly, int lz, const Vec3 &r = Vec3(0, 0, 0));

  GTO(int l = 0, int m = 0, const Vec3 &r = Vec3(0, 0, 0));

  void addPrimitive(double c, double a);

  virtual string toString();
  virtual void przesun(Vec3 &);

  virtual int getL();

  virtual Vec3 getR();

  vector<double> &getC();
  vector<double> &getA();

  inline double x() { return r.x; }
  inline double y() { return r.y; }
  inline double z() { return r.z; }
};

#endif
