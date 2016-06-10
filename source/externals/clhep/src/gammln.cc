// $Id:$
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                         --- HepStat::gammln ---
//                      method implementation file
// -----------------------------------------------------------------------

// =======================================================================
// M. Fischler    - moved the gammln from RandPoisson to here.  01/26/00
// =======================================================================

#include "CLHEP/Random/Stat.h"
#include <cmath>

namespace CLHEP {

double HepStat::gammln(double xx) {

// Returns the value ln(Gamma(xx) for xx > 0.  Full accuracy is obtained for
// xx > 1. For 0 < xx < 1. the reflection formula (6.1.4) can be used first.
// (Adapted from Numerical Recipes in C.  Relative to that routine, this 
// subtracts one from x at the very start, and in exchange does not have to 
// divide ser by x at the end.  The results are formally equal, and practically
// indistinguishable.)

  static const double cof[6] = {76.18009172947146,-86.50532032941677,
                             24.01409824083091, -1.231739572450155,
                             0.1208650973866179e-2, -0.5395239384953e-5};
  int j;
  double x = xx - 1.0;
  double tmp = x + 5.5;
  tmp -= (x + 0.5) * std::log(tmp);
  double ser = 1.000000000190015;

  for ( j = 0; j <= 5; j++ ) {
    x += 1.0;
    ser += cof[j]/x;
  }
  return -tmp + std::log(2.5066282746310005*ser);
}

}  // namespace CLHEP


