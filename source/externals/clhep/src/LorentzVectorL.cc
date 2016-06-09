// -*- C++ -*-
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This is part of the implementation of the HepLorentzVector class:
// Those methods which might, if coded in other modules, force loading 
// of the LorentzRotation.cc code module.
//

#ifdef GNUPRAGMA
#pragma implementation
#endif

#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/LorentzRotation.h"

namespace CLHEP  {

HepLorentzVector &
HepLorentzVector::operator *= (const HepLorentzRotation & m1) {
  return *this = m1.vectorMultiplication(*this);
}

HepLorentzVector &
HepLorentzVector::transform(const HepLorentzRotation & m1){
  return *this = m1.vectorMultiplication(*this);
}

}  // namespace CLHEP
