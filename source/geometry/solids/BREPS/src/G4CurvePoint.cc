// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4CurvePoint.cc,v 1.3 2000-11-08 14:22:09 gcosmo Exp $
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

G4CurvePoint::G4CurvePoint(const G4CurvePoint& cp)
  : c(cp.c), p(cp.p), u(cp.u), notComputed(cp.notComputed)
{
}

G4CurvePoint& G4CurvePoint::operator=(const G4CurvePoint& cp)
{
  if (&cp == this) return *this;
  c = cp.c;
  p = cp.p;
  u = cp.u;
  notComputed = cp.notComputed;
  
  return *this;
}
