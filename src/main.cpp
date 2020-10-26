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

#include <iostream>
#include <basis/basis.h>
#include "job.h"
#include "error.h"

using std::clog;
using std::endl;

int main(int argc, char **argv)
{
  clog << "Hartree-Fock-Roothan, GNU GPLv3 " << endl
       << "Copyright (C) 2009, Sylwester Wysocki <sw143@wp.pl>" << endl
       << "Source code available at https://github.com/dzik143/hfr-cpp" << endl
       << endl;

  if (argc != 2)
  {
    std::cout << "Usage: hfr <file.inp>" << endl;
  }
  else
  {
    // Init basis and geometry classes.
    Basis::checkBasis();
    Geometry::loadSymbols();

    // Create new job object.
    Job job(argv[1]);

    try
    {
      job.run();
    }
    catch (string e)
    {
      std::cerr << e << endl;
    }
    catch (Error e)
    {
      std::cerr << e;
    }
  }
};
