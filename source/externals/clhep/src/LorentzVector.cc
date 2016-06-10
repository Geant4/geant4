// -*- C++ -*-
// $Id:$
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This is the implementation of that portion of the HepLorentzVector class
// which was in the original CLHEP and which does not force loading of either
// Rotation.cc or LorentzRotation.cc
//

#ifdef GNUPRAGMA
#pragma implementation
#endif

#include "CLHEP/Vector/LorentzVector.h"

#include <iostream>

namespace CLHEP  {

double HepLorentzVector::tolerance = 
				Hep3Vector::ToleranceTicks * 2.22045e-16;
double HepLorentzVector::metric = 1.0;

double HepLorentzVector::operator () (int i) const {
  switch(i) {
  case X:
  case Y:
  case Z:
    return pp(i);
  case T:
    return e();
  default:
    std::cerr << "HepLorentzVector subscripting: bad index (" << i << ")"
		 << std::endl;
  }
  return 0.;
}  

double & HepLorentzVector::operator () (int i) {
  static double dummy;
  switch(i) {
  case X:
  case Y:
  case Z:
    return pp(i);
  case T:
    return ee;
  default:
    std::cerr
      << "HepLorentzVector subscripting: bad index (" << i << ")"
      << std::endl;
    return dummy;
  }
}

HepLorentzVector & HepLorentzVector::boost
				(double bx, double by, double bz){
  double b2 = bx*bx + by*by + bz*bz;
  double ggamma = 1.0 / std::sqrt(1.0 - b2);
  double bp = bx*x() + by*y() + bz*z();
  double gamma2 = b2 > 0 ? (ggamma - 1.0)/b2 : 0.0;

  setX(x() + gamma2*bp*bx + ggamma*bx*t());
  setY(y() + gamma2*bp*by + ggamma*by*t());
  setZ(z() + gamma2*bp*bz + ggamma*bz*t());
  setT(ggamma*(t() + bp));
  return *this;
}

HepLorentzVector & HepLorentzVector::rotateX(double a) {
  pp.rotateX(a); 
  return *this; 
}
HepLorentzVector & HepLorentzVector::rotateY(double a) { 
  pp.rotateY(a); 
  return *this; 
}
HepLorentzVector & HepLorentzVector::rotateZ(double a) { 
  pp.rotateZ(a); 
  return *this; 
}

HepLorentzVector & HepLorentzVector::rotateUz(const Hep3Vector &v1) {
  pp.rotateUz(v1);
  return *this;
}

std::ostream & operator<< (std::ostream & os, const HepLorentzVector & v1)
{
  return os << "(" << v1.x() << "," << v1.y() << "," << v1.z()
	    << ";" << v1.t() << ")";
}

std::istream & operator>> (std::istream & is, HepLorentzVector & v1) {

// Required format is ( a, b, c; d ) that is, four numbers, preceded by
// (, followed by ), components of the spatial vector separated by commas,
// time component separated by semicolon. The four numbers are taken
// as x, y, z, t.

  double x, y, z, t;
  char c;

  is >> std::ws >> c;
    // ws is defined to invoke eatwhite(istream & )
    // see (Stroustrup gray book) page 333 and 345.
  if (is.fail() || c != '(' ) {
    std::cerr << "Could not find required opening parenthesis "
	      << "in input of a HepLorentzVector" << std::endl;
    return is;
  }

  is >> x >> std::ws >> c;
  if (is.fail() || c != ',' ) {
    std::cerr << "Could not find x value and required trailing comma "
	      << "in input of a HepLorentzVector" << std::endl; 
    return is;
  }

  is >> y >> std::ws >> c;
  if (is.fail() || c != ',' ) {
    std::cerr << "Could not find y value and required trailing comma "
              <<  "in input of a HepLorentzVector" << std::endl;
    return is;
  }

  is >> z >> std::ws >> c;
  if (is.fail() || c != ';' ) {
    std::cerr << "Could not find z value and required trailing semicolon "
		 <<  "in input of a HepLorentzVector" << std::endl;
    return is;
  }

  is >> t >> std::ws >> c;
  if (is.fail() || c != ')' ) {
    std::cerr << "Could not find t value and required close parenthesis "
		 << "in input of a HepLorentzVector" << std::endl;
    return is;
  }

  v1.setX(x);
  v1.setY(y);
  v1.setZ(z);
  v1.setT(t);
  return is;
}

// The following were added when ZOOM classes were merged in:

HepLorentzVector & HepLorentzVector::operator /= (double c) {
//  if (c == 0) {
//    std::cerr << "HepLorentzVector::operator /=() - "
//      << "Attempt to do LorentzVector /= 0 -- \n"
//      << "division by zero would produce infinite or NAN components"
//      << std::endl;
//  }
  double oneOverC = 1.0/c;
  pp *= oneOverC;
  ee *= oneOverC;
  return *this;
} /* w /= c */

HepLorentzVector operator / (const HepLorentzVector & w, double c) {
//  if (c == 0) {
//    std::cerr << "HepLorentzVector::operator /() - "
//      << "Attempt to do LorentzVector / 0 -- \n"
//      << "division by zero would produce infinite or NAN components"
//      << std::endl;
//  }
  double oneOverC = 1.0/c;
  return HepLorentzVector (w.getV() * oneOverC,
                        w.getT() * oneOverC);
} /* LV = w / c */

Hep3Vector HepLorentzVector::boostVector() const {
  if (ee == 0) {
    if (pp.mag2() == 0) {
      return Hep3Vector(0,0,0);
    } else {
      std::cerr << "HepLorentzVector::boostVector() - "
        << "boostVector computed for LorentzVector with t=0 -- infinite result"
        << std::endl;
      return pp/ee;
    }
  }
  if (restMass2() <= 0) {
    std::cerr << "HepLorentzVector::boostVector() - "
      << "boostVector computed for a non-timelike LorentzVector " << std::endl;
        // result will make analytic sense but is physically meaningless
  }
  return pp * (1./ee);
} /* boostVector */


HepLorentzVector & HepLorentzVector::boostX (double bbeta){
  double b2 = bbeta*bbeta;
  if (b2 >= 1) {
    std::cerr << "HepLorentzVector::boostX() - "
      << "boost along X with beta >= 1 (speed of light) -- \n"
      << "no boost done" << std::endl;
  } else {
    double ggamma = std::sqrt(1./(1-b2));
    double tt = ee;
    ee = ggamma*(ee + bbeta*pp.getX());
    pp.setX(ggamma*(pp.getX() + bbeta*tt));
  }
  return *this;
} /* boostX */

HepLorentzVector & HepLorentzVector::boostY (double bbeta){
  double b2 = bbeta*bbeta;
  if (b2 >= 1) {
    std::cerr << "HepLorentzVector::boostY() - "
      << "boost along Y with beta >= 1 (speed of light) -- \n"
      << "no boost done" << std::endl;
  } else {
    double ggamma = std::sqrt(1./(1-b2));
    double tt = ee;
    ee = ggamma*(ee + bbeta*pp.getY());
    pp.setY(ggamma*(pp.getY() + bbeta*tt));
  }
  return *this;
} /* boostY */

HepLorentzVector & HepLorentzVector::boostZ (double bbeta){
  double b2 = bbeta*bbeta;
  if (b2 >= 1) {
    std::cerr << "HepLorentzVector::boostZ() - "
      << "boost along Z with beta >= 1 (speed of light) -- \n"
      << "no boost done" << std::endl;
  } else {
    double ggamma = std::sqrt(1./(1-b2));
    double tt = ee;
    ee = ggamma*(ee + bbeta*pp.getZ());
    pp.setZ(ggamma*(pp.getZ() + bbeta*tt));
  }
  return *this;
} /* boostZ */

double HepLorentzVector::setTolerance ( double tol ) {
// Set the tolerance for two LorentzVectors to be considered near each other
  double oldTolerance (tolerance);
  tolerance = tol;
  return oldTolerance;
}

double HepLorentzVector::getTolerance ( ) {
// Get the tolerance for two LorentzVectors to be considered near each other
  return tolerance;
}

}  // namespace CLHEP
