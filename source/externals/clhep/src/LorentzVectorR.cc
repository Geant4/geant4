// -*- C++ -*-
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This is part of the implementation of the HepLorentzVector class:
// Those methods which might, if coded in LorentzVector.cc, force loading 
// of the Rotation.cc code module.
//

#ifdef GNUPRAGMA
#pragma implementation
#endif

#include "CLHEP/Vector/LorentzVector.h"

namespace CLHEP  {

HepLorentzVector &  HepLorentzVector::rotate(double a, const Hep3Vector &v1) {
  pp.rotate(a,v1);
  return *this;
}

HepLorentzVector & HepLorentzVector::rotate ( const Hep3Vector & aaxis, 
					      double ddelta )		{
  pp.rotate (aaxis, ddelta);
  return *this;
}

HepLorentzVector & HepLorentzVector::rotate ( const HepAxisAngle & ax ) {
  pp.rotate (ax);
  return *this;
}

HepLorentzVector & HepLorentzVector::rotate ( const HepEulerAngles & e1 ) {
  pp.rotate (e1);
  return *this;
}

HepLorentzVector & HepLorentzVector::rotate ( double phi1,
		                              double theta1,
                		              double psi1 ) {
  pp.rotate (phi1, theta1, psi1);
  return *this;
}

HepLorentzVector rotationOf (const HepLorentzVector & vec, 
			     const Hep3Vector & aaxis,  
			     double ddelta) {
  HepLorentzVector vv (vec);
  return vv.rotate (aaxis, ddelta);
}

HepLorentzVector rotationOf 
	(const HepLorentzVector & vec, const HepAxisAngle &ax ) {
  HepLorentzVector vv (vec);
  return vv.rotate (ax);
}

HepLorentzVector rotationOf
	(const HepLorentzVector & vec, const HepEulerAngles &e1 ) {
  HepLorentzVector vv (vec);
  return vv.rotate (e1);
}

HepLorentzVector rotationOf (const HepLorentzVector & vec, 
				   double phi1,
				   double theta1,
				   double psi1) {
  HepLorentzVector vv (vec);
  return vv.rotate (phi1, theta1, psi1);
}

}  // namespace CLHEP
