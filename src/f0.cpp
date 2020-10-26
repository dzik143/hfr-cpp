#include <cmath>
#include "f0.h"

//=================================================================
//==                                                             ==
//==   CALCULATES THE ERROR FUNCTION ACCORDING TO A RATIONAL     ==
//==   APPROXIMATION FROM M. ABRAMOWITZ AND I.A. STEGUN,         ==
//==   ABSOLUTE ERROR IS LESS THAN 1.5*10**(-7)                  ==
//==   CAN BE REPLACED BY A BUILT-IN FUNCTION ON SOME MACHINES   ==
//==                                                             ==
//=================================================================

double F0(double ARG)
{
  // Ported from FORTRAN 77 code.
  double P = 0.3275911;
  double A[] = {0.254829592, -0.284496736, 1.421413741, -1.453152027, 1.061405429};
  double SARG;
  double ERF;
  double T, TN, POLY;

  // Limit for arg->0.
  if (ARG < 1.0e-6)
  {
    return 1.0 - ARG / 3.0;
  }

  // General case.
  SARG = sqrt(ARG);
  T    = 1.0 / (1.0 + P * SARG);
  TN   = T;
  POLY = A[0] * TN;

  for (int i = 1; i < 5; i++)
  {
    TN   = TN * T;
    POLY = POLY + A[i] * TN;
  }

  ERF = 1.0 - POLY * exp(-ARG);

  return 0.886226925452758014 / SARG * ERF;
}
