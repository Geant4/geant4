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

  const double eps = 1e-15;

  double Ux = rzy - ryz;
  double Uy = rxz - rzx;
  double Uz = ryx - rxy;
  if (std::abs(Ux) < eps && std::abs(Uy) < eps && std::abs(Uz) < eps) {

    double cosdelta = (rxx + ryy + rzz - 1.0) / 2.0;
    if (cosdelta > 0.0) return Hep3Vector(0,0,1); // angle = 0, any axis is good

    double mxx = (rxx + 1)/2;
    double myy = (ryy + 1)/2;
    double mzz = (rzz + 1)/2;
    double mxy = (rxy + ryx)/4;
    double mxz = (rxz + rzx)/4;
    double myz = (ryz + rzy)/4;
    double x, y, z;

    if (mxx > ryy && mxx > rzz) {
      x = std::sqrt(mxx);
      if (rzy - ryz < 0) x = -x;
      y = mxy/x;
      z = mxz/x;
      return  Hep3Vector( x, y, z ).unit();
    } else if (myy > mzz) {
      y = std::sqrt(myy);
      if (rxz - rzx < 0) y = -y;
      x = mxy/y;
      z = myz/y;
      return  Hep3Vector( x, y, z ).unit();
    } else {
      z = std::sqrt(mzz);
      if (ryx - rxy < 0) z = -z;
      x = mxz/z;
      y = myz/z;
      return  Hep3Vector( x, y, z ).unit();
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
