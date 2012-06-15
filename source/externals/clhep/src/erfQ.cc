// $Id:$
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                          --- erfQ ---
//                      method implementation file
// -----------------------------------------------------------------------

// Contains methods that do not depend on large tables.
//
// erfQ		  (double x)

// =======================================================================
// M Fischler	  - Created 1/26/00.
//
// =======================================================================

#include "CLHEP/Random/Stat.h"
#include <cmath>

namespace CLHEP {

double HepStat::erfQ (double x) {
//
// erfQ is accurate to 7 places.
// See Numerical Recipes P 221
//

  double t, z, erfc; 

  z = std::abs(x);
  t = 1.0/(1.0+.5*z);

  erfc= t*std::exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
	t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
	t*(-0.82215223+t*0.17087277 ))) ))) )));

  // (The derivation of this formula should be obvious.)

  if ( x < 0 ) erfc = 2.0 - erfc;

  return 1 - erfc;

}


}  // namespace CLHEP

