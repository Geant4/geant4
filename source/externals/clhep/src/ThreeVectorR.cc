// -*- C++ -*-
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This is the implementation of those methods of the Hep3Vector class which
// require linking of the HepRotation class.  These methods have been broken 
// out of ThreeVector.cc.
//

#ifdef GNUPRAGMA
#pragma implementation
#endif

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

namespace CLHEP  {

Hep3Vector & Hep3Vector::operator *= (const HepRotation & m1) {
  return *this = m1 * (*this);
}

Hep3Vector & Hep3Vector::transform(const HepRotation & m1) {
  return *this = m1 * (*this);
}

Hep3Vector & Hep3Vector::rotate(double angle1, const Hep3Vector & aaxis){
  HepRotation trans;
  trans.rotate(angle1, aaxis);
  operator*=(trans);
  return *this;
}

}  // namespace CLHEP
