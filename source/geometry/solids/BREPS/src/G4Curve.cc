// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Curve.cc,v 1.3 2000-08-28 08:57:57 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4Curve.cc
//
// ----------------------------------------------------------------------

#include "G4Curve.hh"

G4Curve::G4Curve()
 : bounded(false), bBox(G4BoundingBox3D::space), sameSense(true)
{
}

G4Curve::~G4Curve()
{
}

G4String G4Curve::GetEntityType() const
{
  return G4String("G4Curve");
}

const char* G4Curve::Name() const
{
  return "G4Curve";
}

void G4Curve::SetParentSrfPtr(const G4Surface*)
{
}
