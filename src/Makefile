################################################################################
#									       #
# Copyright (C) 2009 Sylwester Wysocki <sw143@wp.pl>			       #
#									       #
# Source code available at: https://github.com/dzik143/hfr-cpp		       #
#									       #
# This program is free software: you can redistribute it and/or modify	       #
# it under the terms of the GNU General Public License as published by	       #
# the Free Software Foundation, either version 3 of the License, or	       #
# (at your option) any later version.					       #
#									       #
# This program is distributed in the hope that it will be useful,	       #
# but WITHOUT ANY WARRANTY; without even the implied warranty of	       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the 	       #
# GNU General Public License for more details.				       #
#									       #
# You should have received a copy of the GNU General Public License	       #
# along with this program. If not, see <http://www.gnu.org/licenses/>	       #
#									       #
################################################################################

CXXSRC = \
  f0.cpp \
  error.cpp \
  job.cpp \
  main.cpp \
  algebra/vector.cpp \
  algebra/matrix.cpp \
  algebra/tensor4.cpp \
  basis/basis.cpp \
  basis/gto.cpp \
  geometry/vec3.cpp \
  geometry/atom.cpp \
  geometry/geometry.cpp \
  methods/hfr.cpp \
  psi/psi.cpp

CXXOBJ	 = $(CXXSRC:.cpp=.o)
CFLAGS	 = -O2
CXX	 = g++
CXXFLAGS = -I.

all: hfr

$(CXXOBJ) : %.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< $(CFLAGS) -o $@

hfr: $(CXXOBJ)
	$(CXX) $(CXXOBJ) -o $@

clean:
	rm -f $(CXXOBJ)
	rm -f *.log
	rm -f *.a
	rm -f *.def
	rm -f hfr
	rm -f hfr.exe
