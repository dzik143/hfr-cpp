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

#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>

class Matrix
{
  friend std::istream &operator>>(std::istream &, Matrix &);
  friend std::ostream &operator<<(std::ostream &, Matrix &);

protected:

  double *data;

  int rowsCnt;
  int colsCnt;

public:

  Matrix(int rowsCnt = 0, int colsCnt = 0);
  Matrix(int rowsCnt, int colsCnt, const double **initialData);

  Matrix(const Matrix &);

  ~Matrix();

  void loadIdentity();

  double &operator()(int, int);      // operatory dostepu
  double operator()(int, int) const; // do komorki (i,j)

  double &operator[](int i);      // operatory dostepu
  double operator[](int i) const; // do 1-szej kolumny

  Matrix &operator=(const Matrix &);
  Matrix operator+(const Matrix &);
  Matrix operator-(const Matrix &);
  Matrix operator*(const Matrix &);

  Matrix operator*(const double);

  Matrix transpone(); // transponowanie
  Matrix invert();    // odwraca macierz
  void power(double); // potegowanie macierzy

  int szer();   // zwraca szerokosc
  int wys();    // zwraca wysokosc
  bool empty(); // true, gdy pusta

  void diagonalize(Matrix &, Matrix &); // diagnonalizacja
  double sgn(double);                   // zwraca znak argumentu

  Matrix pierw(); // pierwiastek
};

std::istream &operator>>(std::istream &, Matrix &); // wczytywanie ze strumienia
std::ostream &operator<<(std::ostream &, Matrix &); // wypisywanie na strumienia

#endif
