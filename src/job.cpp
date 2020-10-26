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

#include "job.h"

using std::clog;
using std::endl;

Job::Job()
{
  this->baza = NULL;
  this->psi  = NULL;
};

Job::~Job()
{
  if (this->baza != NULL)
  {
    delete this->baza;
  }

  if (this->psi != NULL)
  {
    delete this->psi;
  }

  if (this->input)
  {
    input->close();
    delete this->input;
  }
};

Job::Job(const string &fname)
{
  this->baza  = NULL;
  this->psi   = NULL;
  this->input = new std::ifstream(fname);

  clog << "Using job file '" << fname << "'" << endl;
}

void Job::run()
{
  while (!(*input).eof())
  {
    // Read next command.
    std::string komenda;
    *input >> std::ws >> komenda;

    // ========================
    // = Geometry declaration =
    // ========================

    if (komenda == "GEOM")
    {
      geom.loadGeom(*input);
    }

    // =========================
    // = Basis set declaration =
    // =========================

    else if (komenda == "BASIS")
    {
      string bazaName;
      string bazaType;
      *input >> bazaName;

      bazaType = Basis::basisType[bazaName];

      clog << "Using basis set " << bazaName << "..." << endl;

      if (bazaType == "GTO")
      {
        // Create new GTO basis set for our system.
        baza = new BasisGTO();

        // Load basis set.
        baza->loadBasis(Basis::basisFile[bazaName]);

        // Create atomic orbitals for given geometry (molecule)
        // using basis functions.
        baza->initAO(geom);
      }
      else
      {
        throw string("invalid basis set: ") + bazaName;
      }
    }

    // ===========================
    // = Run Hartree-Fock method =
    // ===========================

    else if (komenda == "HFR")
    {
      HartreeFockRoothan hfr;
      hfr.solve(*baza, geom);
    }

    // ===========================
    // = Error - unknown command =
    // ===========================

    else
    {
      throw string("invalid command: ") + komenda;
    }
  }
}
