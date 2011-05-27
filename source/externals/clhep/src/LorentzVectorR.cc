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

HepLorentzVector &  HepLorentzVector::rotate(double a, const Hep3Vector &v) {
  pp.rotate(a,v);
  return *this;
}

HepLorentzVector & HepLorentzVector::rotate ( const Hep3Vector & axis, 
					      double delta )		{
  pp.rotate (axis, delta);
  return *this;
}

HepLorentzVector & HepLorentzVector::rotate ( const HepAxisAngle & ax ) {
  pp.rotate (ax);
  return *this;
}

HepLorentzVector & HepLorentzVector::rotate ( const HepEulerAngles & e ) {
  pp.rotate (e);
  return *this;
}

HepLorentzVector & HepLorentzVector::rotate ( double phi,
		                              double theta,
                		              double psi ) {
  pp.rotate (phi, theta, psi);
  return *this;
}

HepLorentzVector rotationOf (const HepLorentzVector & vec, 
			     const Hep3Vector & axis,  
			     double delta) {
  HepLorentzVector vv (vec);
  return vv.rotate (axis, delta);
}

HepLorentzVector rotationOf 
	(const HepLorentzVector & vec, const HepAxisAngle &ax ) {
  HepLorentzVector vv (vec);
  return vv.rotate (ax);
}

HepLorentzVector rotationOf
	(const HepLorentzVector & vec, const HepEulerAngles &e ) {
  HepLorentzVector vv (vec);
  return vv.rotate (e);
}

HepLorentzVector rotationOf (const HepLorentzVector & vec, 
				   double phi,
				   double theta,
				   double psi) {
  HepLorentzVector vv (vec);
  return vv.rotate (phi, theta, psi);
}

}  // namespace CLHEP
