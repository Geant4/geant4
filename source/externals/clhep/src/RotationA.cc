// -*- C++ -*-
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This is the implementation of those methods of the HepRotation class which
// were introduced when ZOOM PhysicsVectors was merged in, and which involve 
// the angle/axis representation of a Rotation.
//

#ifdef GNUPRAGMA
#pragma implementation
#endif

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include <iostream>
#include <cmath>

namespace CLHEP  {

// ----------  Constructors and Assignment:

// axis and angle

HepRotation & HepRotation::set( const Hep3Vector & aaxis, double ddelta ) {

  double sinDelta = std::sin(ddelta), cosDelta = std::cos(ddelta);
  double oneMinusCosDelta = 1.0 - cosDelta;

  Hep3Vector u = aaxis.unit();

  double uX = u.getX();
  double uY = u.getY();
  double uZ = u.getZ();

  rxx = oneMinusCosDelta * uX * uX  +  cosDelta;
  rxy = oneMinusCosDelta * uX * uY  -  sinDelta * uZ;
  rxz = oneMinusCosDelta * uX * uZ  +  sinDelta * uY;

  ryx = oneMinusCosDelta * uY * uX  +  sinDelta * uZ;
  ryy = oneMinusCosDelta * uY * uY  +  cosDelta;
  ryz = oneMinusCosDelta * uY * uZ  -  sinDelta * uX;

  rzx = oneMinusCosDelta * uZ * uX  -  sinDelta * uY;
  rzy = oneMinusCosDelta * uZ * uY  +  sinDelta * uX;
  rzz = oneMinusCosDelta * uZ * uZ  +  cosDelta;

  return  *this;

} // HepRotation::set(axis, delta)

HepRotation::HepRotation ( const Hep3Vector & aaxis, double ddelta ) 
{
  set( aaxis, ddelta );
}  
HepRotation & HepRotation::set( const HepAxisAngle & ax ) {
  return  set ( ax.axis(), ax.delta() );
}
HepRotation::HepRotation ( const HepAxisAngle & ax ) 
{
  set ( ax.axis(), ax.delta() );
}

double    HepRotation::delta() const {

  double cosdelta = (rxx + ryy + rzz - 1.0) / 2.0;
  if (cosdelta > 1.0) {
    return 0;
  } else if (cosdelta < -1.0) {
    return CLHEP::pi;
  } else {
    return  std::acos( cosdelta ); // Already safe due to the cosdelta > 1 check
  }

} // delta()

Hep3Vector HepRotation::axis () const {

  // Determine 2*std::sin(delta) times the u components (I call this uX, uY, Uz)
  // Normalization is not needed; it will be done when returning the 3-Vector

  double  Uz = ryx - rxy;
  double  Uy = rxz - rzx;
  double  Ux = rzy - ryz;

  if ( (Uz==0) && (Uy==0) && (Ux==0) ) {
    if        ( rzz>0 ) {
      return Hep3Vector(0,0,1);
    } else if ( ryy>0 ) {
      return Hep3Vector(0,1,0);
    } else {
      return Hep3Vector(1,0,0);
    }
  } else {
    return  Hep3Vector( Ux, Uy, Uz ).unit();
  }

} // axis()

HepAxisAngle HepRotation::axisAngle() const {

  return HepAxisAngle (axis(), delta());

} // axisAngle() 


void HepRotation::setAxis (const Hep3Vector & aaxis) {
  set ( aaxis, delta() );
}

void HepRotation::setDelta (double ddelta) {
  set ( axis(), ddelta );
}

}  // namespace CLHEP
