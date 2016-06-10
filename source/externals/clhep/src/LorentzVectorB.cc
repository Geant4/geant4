// -*- C++ -*-
// $Id:$
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This is the implementation of the HepLorentzVector class:
// Those methods originating in ZOOM dealing with simple boosts and rotations.
// Use of one of these methods will not force loading of the HepRotation or
// HepLorentzRotation class.
//

#ifdef GNUPRAGMA
#pragma implementation
#endif

#include "CLHEP/Vector/LorentzVector.h"

namespace CLHEP  {

//-*********
// rotationOf()
//-*********

// Each of these is a shell over a rotate method.

HepLorentzVector rotationXOf
	(const HepLorentzVector & vec, double phi){
  HepLorentzVector vv (vec);
  return vv.rotateX (phi);
}

HepLorentzVector rotationYOf
	(const HepLorentzVector & vec, double phi){
  HepLorentzVector vv (vec);
  return vv.rotateY (phi);
}

HepLorentzVector rotationZOf
	(const HepLorentzVector & vec, double phi){
  HepLorentzVector vv (vec);
  return vv.rotateZ (phi);
}

//-********
// boost
//-********

HepLorentzVector & HepLorentzVector::boost 
			( const Hep3Vector & aaxis,  double bbeta ) {
  if (bbeta==0) {
    return *this; // do nothing for a 0 boost
  }
  double r2 = aaxis.mag2();
  if ( r2 == 0 ) {
    std::cerr << "HepLorentzVector::boost() - "
      << "A zero vector used as axis defining a boost -- no boost done"
      << std::endl;
    return *this;
  } 
  double b2 = bbeta*bbeta;
  if (b2 >= 1) {
    std::cerr << "HepLorentzVector::boost() - "
      << "LorentzVector boosted with beta >= 1 (speed of light) -- \n"
      << "no boost done" << std::endl;
  } else {
    Hep3Vector u = aaxis.unit();
    double ggamma = std::sqrt(1./(1.-b2));
    double betaDotV = u.dot(pp)*bbeta;
    double tt = ee;

    ee = ggamma * (tt + betaDotV);
    pp += ( ((ggamma-1)/b2)*betaDotV*bbeta + ggamma*bbeta*tt ) * u;
    // Note:  I have verified the behavior of this even when beta is very
    //        small -- (gamma-1)/b2 becomes inaccurate by O(1), but it is then
    //        multiplied by O(beta**2) and added to an O(beta) term, so the
    //        inaccuracy does not affect the final result.
  }
  return *this;
} /* boost (axis, beta) */

}  // namespace CLHEP
