// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4CurvePoint.cc,v 1.2 2000-08-28 08:57:57 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4CurvePoint.cc
//
// ----------------------------------------------------------------------

#include "G4CurvePoint.hh"

const G4int G4CurvePoint::pFlag= 1;
const G4int G4CurvePoint::uFlag= 2;
const G4int G4CurvePoint::allFlags= 0xFF; // lots of bits...

G4CurvePoint::G4CurvePoint()
 : c(0)
{
}

G4CurvePoint::G4CurvePoint(G4Curve& c0)
{
  Init(c0);
}

G4CurvePoint::~G4CurvePoint()
{
}
