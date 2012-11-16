// -*- C++ -*-
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This is the implementation of methods of the HepRotationX class which
// were introduced when ZOOM PhysicsVectors was merged in.
//

#ifdef GNUPRAGMA
#pragma implementation
#endif

#include "CLHEP/Vector/RotationX.h"
#include "CLHEP/Vector/AxisAngle.h"
#include "CLHEP/Vector/EulerAngles.h"
#include "CLHEP/Vector/LorentzRotation.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include <cmath>
#include <stdlib.h>
#include <iostream>

namespace CLHEP  {

static inline double safe_acos (double x) {
  if (std::abs(x) <= 1.0) return std::acos(x);
  return ( (x>0) ? 0 : CLHEP::pi );
}

HepRotationX::HepRotationX(double ddelta) : 
		its_d(proper(ddelta)), its_s(std::sin(ddelta)), its_c(std::cos(ddelta))
{}

HepRotationX & HepRotationX::set ( double ddelta ) {
  its_d = proper(ddelta);
  its_s = std::sin(its_d);
  its_c = std::cos(its_d);
  return *this;
}

double  HepRotationX::phi() const {
  if ( (its_d > 0) && (its_d < CLHEP::pi) ) {
    return CLHEP::pi;
  } else {
    return 0.0;
  }
}  // HepRotationX::phi()

double  HepRotationX::theta() const {
  return  std::fabs( its_d );
}  // HepRotationX::theta()

double  HepRotationX::psi() const {
  if ( (its_d > 0) && (its_d < CLHEP::pi) ) {
    return CLHEP::pi;
  } else {
    return 0.0;
  }
}  // HepRotationX::psi()

HepEulerAngles HepRotationX::eulerAngles() const {
  return HepEulerAngles(  phi(), theta(),  psi() );
}  // HepRotationX::eulerAngles()


// From the defining code in the implementation of CLHEP (in Rotation.cc)
// it is clear that thetaX, phiX form the polar angles in the original
// coordinate system of the new X axis (and similarly for phiY and phiZ).
//
// This code is taken directly from the original CLHEP. However, there are as
// shown opportunities for significant speed improvement.

double HepRotationX::phiX() const {
  return (yx() == 0.0 && xx() == 0.0) ? 0.0 : std::atan2(yx(),xx());
  		// or ---- return 0;
}

double HepRotationX::phiY() const {
  return (yy() == 0.0 && xy() == 0.0) ? 0.0 : std::atan2(yy(),xy());
		// or ----  return (yy() == 0.0) ? 0.0 : std::atan2(yy(),xy());
}

double HepRotationX::phiZ() const {
  return (yz() == 0.0 && xz() == 0.0) ? 0.0 : std::atan2(yz(),xz());
		// or ----  return (yz() == 0.0) ? 0.0 : std::atan2(yz(),xz());
}

double HepRotationX::thetaX() const {
  return safe_acos(zx());
		// or ----  return CLHEP::halfpi;
}

double HepRotationX::thetaY() const {
  return safe_acos(zy());
}

double HepRotationX::thetaZ() const {
  return safe_acos(zz());  
		// or ---- return d;
}

void HepRotationX::setDelta ( double ddelta ) {
  set(ddelta);
}

void HepRotationX::decompose
	(HepAxisAngle & rotation, Hep3Vector & boost) const {
  boost.set(0,0,0);
  rotation = axisAngle();
}

void HepRotationX::decompose
	(Hep3Vector & boost, HepAxisAngle & rotation) const {
  boost.set(0,0,0);
  rotation = axisAngle();
}

void HepRotationX::decompose
        (HepRotation & rotation, HepBoost & boost) const {
  boost.set(0,0,0);
  rotation = HepRotation(*this);
} 

void HepRotationX::decompose
        (HepBoost & boost, HepRotation & rotation) const {
  boost.set(0,0,0);
  rotation = HepRotation(*this);
}

double HepRotationX::distance2( const HepRotationX & r  ) const {
  double answer = 2.0 * ( 1.0 - ( its_s * r.its_s + its_c * r.its_c ) ) ;
  return (answer >= 0) ? answer : 0;
}

double HepRotationX::distance2( const HepRotation & r  ) const {
  double sum =        r.xx() + 
                		  yy() * r.yy() + yz() * r.yz()
                		+ zy() * r.zy() + zz() * r.zz();
  double answer = 3.0 - sum;
  return (answer >= 0 ) ? answer : 0;
}

double HepRotationX::distance2( const HepLorentzRotation & lt  ) const {
  HepAxisAngle a; 
  Hep3Vector   b;
  lt.decompose(b, a);
  double bet = b.beta();
  double bet2 = bet*bet;
  HepRotation r(a);
  return bet2/(1-bet2) + distance2(r);
}

double HepRotationX::distance2( const HepBoost & lt ) const {
  return distance2( HepLorentzRotation(lt));
}

double HepRotationX::howNear( const HepRotationX & r ) const {
  return std::sqrt(distance2(r));
}
double HepRotationX::howNear( const HepRotation & r ) const {
  return std::sqrt(distance2(r));
}
double HepRotationX::howNear( const HepBoost & b ) const {
  return std::sqrt(distance2(b));
}
double HepRotationX::howNear( const HepLorentzRotation & lt ) const {
  return std::sqrt(distance2(lt));
}
bool HepRotationX::isNear(const HepRotationX & r,double epsilon)const{
  return (distance2(r) <= epsilon*epsilon);
}
bool HepRotationX::isNear(const HepRotation & r,double epsilon) const{
  return (distance2(r) <= epsilon*epsilon);
}
bool HepRotationX::isNear( const HepBoost & lt,double epsilon) const {
  return (distance2(lt) <= epsilon*epsilon);
}

bool HepRotationX::isNear( const HepLorentzRotation & lt,
                                     double epsilon ) const {
  return (distance2(lt) <= epsilon*epsilon);
}

double HepRotationX::norm2() const {
  return 2.0 - 2.0 * its_c;
}

std::ostream & HepRotationX::print( std::ostream & os ) const {
  os << "\nRotation about X (" << its_d << 
		") [cos d = " << its_c << " sin d = " << its_s << "]\n";
  return os;
}

}  // namespace CLHEP

