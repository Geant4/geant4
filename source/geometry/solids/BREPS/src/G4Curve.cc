// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Curve.cc,v 1.5 2000-11-20 17:54:39 gcosmo Exp $
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
 : bBox(G4BoundingBox3D::space), bounded(false), sameSense(true)
{
}

G4Curve::~G4Curve()
{
}

G4Curve::G4Curve(const G4Curve& c)
 : start(c.start), end(c.end), pStart(c.pStart), pEnd(c.pEnd),
   pRange(c.pRange), bounded(c.bounded), sameSense(c.sameSense)
{
}

G4Curve& G4Curve::operator=(const G4Curve& c)
{
  if (&c == this) return *this;
  start     = c.start;
  end       = c.end;
  pStart    = c.pStart;
  pEnd      = c.pEnd;
  pRange    = c.pRange;
  bounded   = c.bounded;
  sameSense = c.sameSense;

  return *this;
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
