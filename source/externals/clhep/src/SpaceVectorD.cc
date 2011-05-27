// -*- C++ -*-
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This is the implementation of the subset of those methods of the Hep3Vector 
// class which originated from the ZOOM SpaceVector class *and* which involve
// the esoteric concepts of polar/azimuthal angular decomposition.
//

#ifdef GNUPRAGMA
#pragma implementation
#endif

#include "CLHEP/Vector/ThreeVector.h"

#include <cmath>

namespace CLHEP  {

//-*********************************************
//			- 6 -
// Decomposition of an angle between two vectors
//
//-*********************************************


double Hep3Vector::polarAngle (const Hep3Vector & v2) const {
  return std::fabs(v2.getTheta() - getTheta());
} /* polarAngle */

double Hep3Vector::polarAngle (const Hep3Vector & v2,
				const Hep3Vector & ref) const {
  return std::fabs( v2.angle(ref) - angle(ref) );
} /* polarAngle (v2, ref) */

// double Hep3Vector::azimAngle (const Hep3Vector & v2) const 
// is now in the .icc file as deltaPhi(v2)

double Hep3Vector::azimAngle  (const Hep3Vector & v2,
				const Hep3Vector & ref) const {

  Hep3Vector vperp ( perpPart(ref) );
  if ( vperp.mag2() == 0 ) {
    std::cerr << "Hep3Vector::azimAngle() - "
      << "Cannot find azimuthal angle with reference direction parallel to "
      << "vector 1 -- will return zero" << std::endl;
   return 0;
  }

  Hep3Vector v2perp ( v2.perpPart(ref) );
  if ( v2perp.mag2() == 0 ) {
    std::cerr << "Hep3Vector::azimAngle() - "
      << "Cannot find azimuthal angle with reference direction parallel to "
      << "vector 2 -- will return zero" << std::endl;
   return 0;
  }

  double ang = vperp.angle(v2perp);

  // Now compute the sign of the answer:  that of U*(VxV2) or 
  // the equivalent expression V*(V2xU).

  if  ( dot(v2.cross(ref)) >= 0 ) {
    return ang;
  } else {
    return -ang;
  }

	//-| Note that if V*(V2xU) is zero, we want to return 0 or PI
	//-| depending on whether vperp is aligned or antialigned with v2perp.
	//-| The computed angle() expression does this properly.

} /* azimAngle (v2, ref) */

}  // namespace CLHEP
