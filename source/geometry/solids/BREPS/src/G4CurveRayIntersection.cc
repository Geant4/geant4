// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4CurveRayIntersection.cc,v 1.2 2000-08-28 08:57:57 gcosmo Exp $
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
{
  d= kInfinity;
}

G4CurveRayIntersection::G4CurveRayIntersection(G4Curve& c0, const G4Ray& r0)
{
  Init(c0, r0);
}

G4CurveRayIntersection::~G4CurveRayIntersection()
{
}
