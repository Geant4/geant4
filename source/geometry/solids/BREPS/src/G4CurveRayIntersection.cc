// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4CurveRayIntersection.cc,v 1.3 2000-11-08 14:22:10 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4CurveRayIntersection.cc
//
// ----------------------------------------------------------------------

#include "G4CurveRayIntersection.hh"

const G4int G4CurveRayIntersection::dFlag= 4;

G4CurveRayIntersection::G4CurveRayIntersection()
  : r(0), d(kInfinity)
{
}

G4CurveRayIntersection::G4CurveRayIntersection(G4Curve& c0, const G4Ray& r0)
{
  Init(c0, r0);
}

G4CurveRayIntersection::~G4CurveRayIntersection()
{
}

G4CurveRayIntersection::G4CurveRayIntersection(const G4CurveRayIntersection& cr)
  : r(cr.r), d(cr.d)
{
  c = cr.c;
  p = cr.p;
  u = cr.u;
  notComputed = cr.notComputed;
}

G4CurveRayIntersection&
G4CurveRayIntersection::operator=(const G4CurveRayIntersection& cr)
{
  if (&cr == this) return *this;
  c = cr.c;
  p = cr.p;
  u = cr.u;
  r = cr.r;
  d = cr.d;
  notComputed = cr.notComputed;

  return *this;
}
