// -*- C++ -*-
// $Id:$
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This is the implementation of those few parts of the Hep4RotationInterface
// and Hep3RotationInterface classes which are neither inline nor pure virtual.
//

#ifdef GNUPRAGMA
#pragma implementation
#endif

#include "CLHEP/Vector/RotationInterfaces.h"

namespace CLHEP  {

//-******************************
//
// Hep4RotationInterface
//
//-******************************

double Hep4RotationInterface::tolerance = 
			Hep4RotationInterface::ToleranceTicks * 1.0e-08;

double Hep4RotationInterface::getTolerance() {return tolerance;} 
double Hep4RotationInterface::setTolerance( double tol ) {
  double t = tolerance; tolerance = tol; return t;
}


//-******************************
//
// Hep3RotationInterface
//
//-******************************

}  // namespace CLHEP
