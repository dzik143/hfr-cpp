# [ARCHIVE/2008] Hartree-Fock-Roothan (C++)
  - **NAIVE** implementation of [Hartree-Fock method](https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method) in [GTO base](https://en.wikipedia.org/wiki/Gaussian_orbital),
  - **ARCHIVAL** code written in **2009**,
  - port of older [fortran code](https://github.com/dzik143/hfr-f77),
  - **GPLv3 LICENSE** - use for any purpose as long as you keep copyright notice and keep open source.

# Features:
  - **NAIVE SCF** loop - O(N^4) complexity,
  - calculates **ENERGY** and **ORBITALS** (represented by coefficients in the matrix),
  - read molecule **GEOMETRY** from file,
  - **GTO BASIS** with parameters read from file.

# Limitations:
  - Only integrals for [S orbitals](https://en.wikipedia.org/wiki/Gaussian_orbital) are implemented,
  - the number of total electrons in the system **MUST BE EVEN** (see [restricted Hartree-Fock](https://en.wikipedia.org/wiki/Restricted_open-shell_Hartree%E2%80%93Fock)).

# Build:
  ```
  cd src
  make
  ```

# Usage:
  - Create input file containing:
    - molecule geometry in cartessian system,
    - basis name to use (basis must be declared within *basis.lib* file),
    - example *h2.inp* file is presented below:
  ```
GEOM
 H 0.0 0.0 0.0
 H 0.0 0.0 0.742724
END_GEOM

BASIS 6-31G
HFR
```

  - Run HFR code:
  ```
  $./hfr h2.inp
Hartree-Fock-Roothan, GNU GPLv3
Copyright (C) 2009, Sylwester Wysocki <sw143@wp.pl>
Source code available at https://github.com/dzik143/hfr-cpp

Loading basis list.................O.K!
Loading element symbols............O.K!
Using job file 'h2.inp'
Loading geometry...................O.K!

 INPUT GEOMETRY [au]
-------------------------------
  Atom      X      Y      Z
   H   0.00000 0.00000 0.00000
   H   0.00000 0.00000 1.40348
-------------------------------

Using basis set 6-31G...
Loading basis set..................O.K!
Generating atomic orbitals.........O.K!

----------------------------------------------------
   Number of contracted GTO functions: 4
   Number of primitive GTO functions:  8
   Number of one-electron integrals:   64
   Number of two-electron integrals:   1024
   Memory needed to store integrals:   0.004 MB
----------------------------------------------------

Calculating molecular integrals....O.K!
Falling into SCF loop...

iter=1   Etotal=-1.07423828000
iter=2   Etotal=-1.12542019000 delta=5.11819070e-002
iter=3   Etotal=-1.12668272000 delta=1.26253187e-003
iter=4   Etotal=-1.12671112000 delta=2.84026845e-005
iter=5   Etotal=-1.12671175000 delta=6.30249303e-007
iter=6   Etotal=-1.12671177000 delta=1.39565290e-008

SCF convergented.

Molecular orbitals:
   1.120   0.123   0.768   0.327
  -1.347   1.708  -0.686   0.272
  -1.120  -0.123   0.768   0.327
   1.347  -1.708  -0.686   0.272

Orbital energies:
   1.401
   0.238
   0.776
  -0.595
```
