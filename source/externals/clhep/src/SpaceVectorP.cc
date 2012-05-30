// -*- C++ -*-
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// SpaceVector
//
// This is the implementation of the subset of those methods of the Hep3Vector 
// class which originated from the ZOOM SpaceVector class *and* which involve
// intrinsic properties or propeties relative to a second vector.
//

#ifdef GNUPRAGMA
#pragma implementation
#endif

#include "CLHEP/Vector/ThreeVector.h"

#include <cmath>

namespace CLHEP  {

//-********************************
//		- 5 -
// Intrinsic properties of a vector
// and properties relative to a direction
//
//-********************************

double Hep3Vector::beta() const {
  double b = std::sqrt(mag2());
//  if (b >= 1) {
//    std::cerr << "Hep3Vector::beta() - "
//      << "Beta taken for Hep3Vector of at least unit length" << std::endl;
//  }
  return b;
}

double Hep3Vector::gamma() const {
  double bbeta = std::sqrt(mag2());
//  if (bbeta == 1) {
//    std::cerr << "Hep3Vector::gamma() - "
//      << "Gamma taken for Hep3Vector of unit magnitude -- infinite result"
//      << std::endl;
//  }
//  if (bbeta > 1) {
//    std::cerr << "Hep3Vector::gamma() - "
//      << "Gamma taken for Hep3Vector of more than unit magnitude -- \n"
//      << "the sqrt function would return NAN" << std::endl;
//  }
  return 1/std::sqrt(1-bbeta*bbeta);
}

double Hep3Vector::rapidity() const {
//  if (std::fabs(dz) == 1) {
//    std::cerr << "Hep3Vector::rapidity() - "
//      << "Rapidity in Z direction taken for Hep3Vector with |Z| = 1 -- \n"
//      << "the log should return infinity" <, std::endl;
//  }
//  if (std::fabs(dz) > 1) {
//    std::cerr << "Hep3Vector::rapidity() - "
//      << "Rapidity in Z direction taken for Hep3Vector with |Z| > 1 -- \n"
//      << "the log would return a NAN" << std::endl;
//  }
  // Want inverse std::tanh(dz):
  return (.5 * std::log((1+dz)/(1-dz)) );
}

double Hep3Vector::coLinearRapidity() const {
  double b = beta();
//  if (b == 1) {
//    std::cerr << "Hep3Vector::coLinearRapidity() - "
//      << "Co-linear Rapidity taken for Hep3Vector of unit length -- \n"
//      << "the log should return infinity" << std::endl;
//  }
//  if (b > 1) {
//    std::cerr << "Hep3Vector::coLinearRapidity() - "
//      << "Co-linear Rapidity taken for Hep3Vector of more than unit length -- \n"
//      << "the log would return a NAN" << std::endl;
//  }
  // Want inverse std::tanh(b):
  return (.5 * std::log((1+b)/(1-b)) );
}

//-***********************************************
// Other properties relative to a reference vector
//-***********************************************

Hep3Vector Hep3Vector::project (const Hep3Vector & v2) const {
  double mag2v2 = v2.mag2();
  if (mag2v2 == 0) {
    std::cerr << "Hep3Vector::project() - "
      << "Attempt to take projection of vector against zero reference vector"
      << std::endl;
    return project();
  }
  return ( v2 * (dot(v2)/mag2v2) );
}

double Hep3Vector::rapidity(const Hep3Vector & v2) const {
  double vmag = v2.mag();
  if ( vmag == 0 ) {
    std::cerr << "Hep3Vector::rapidity() - "
      << "Rapidity taken with respect to zero vector" << std::endl;
    return 0;    
  }
  double z1 = dot(v2)/vmag;
//  if (std::fabs(z1) >= 1) {
//    std::cerr << "Hep3Vector::rapidity() - "
//      << "Rapidity taken for too large a Hep3Vector "
//      << "-- would return infinity or NAN" << std::endl;
//  }
  // Want inverse std::tanh(z):
  return (.5 * std::log((1+z1)/(1-z1)) );
}

double Hep3Vector::eta(const Hep3Vector & v2) const {
  // Defined as    -std::log ( std::tan ( .5* theta(u) ) );
  //
  // Quicker is to use cosTheta:
  // std::tan (theta/2) = std::sin(theta)/(1 + std::cos(theta))

  double r1   = getR();
  double v2r = v2.mag();
  if ( (r1 == 0) || (v2r == 0) ) {
    std::cerr << "Hep3Vector::eta() - "
      << "Cannot find pseudorapidity of a zero vector relative to a vector"
      << std::endl;
    return 0.;
  }
  double c  = dot(v2)/(r1*v2r);
  if ( c >= 1 ) {
    c = 1; 	//-| We don't want to return NAN because of roundoff
    std::cerr << "Hep3Vector::eta() - "
      << "Pseudorapidity of vector relative to parallel vector -- \n"
      << "will give infinite result" << std::endl;
   			    // We can just go on; tangent will be 0, so
			    // std::log (tangent) will be -INFINITY, so result
			    // will be +INFINITY.
  }
  if ( c <= -1 ) {
    std::cerr << "Hep3Vector::eta() - "
      << "Pseudorapidity of vector relative to anti-parallel vector -- \n"
      << "will give negative infinite result"<< std::endl;
    			//-| We don't want to return NAN because of roundoff
    return ( negativeInfinity() );
			    //  If we just went on, the tangent would be NAN
			    //  so return would be NAN.  But the proper limit
			    // of tan is +Infinity, so the return should be
			    // -INFINITY.
  }

  double tangent = std::sqrt (1-c*c) / ( 1 + c );
  return (- std::log (tangent));

} /* eta (u) */


}  // namespace CLHEP
