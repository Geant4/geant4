// Persistent class describing solid placements for boolean operations
//
// History:
// 10.11.99 Y.Morita, Initial creation

#include "G4PDisplacedSolid.hh"
#include "G4DisplacedSolid.hh"
#include "G4AffineTransform.hh"

#include "G4VSolid.hh"

G4PDisplacedSolid::G4PDisplacedSolid
                       ( HepRef(G4PVSolid) persCostituentSolid,
                         HepRef(G4PAffineTransform) pDirectTransform )
{
  fPtrSolid = persCostituentSolid;
  fDirectTransform = pDirectTransform;
}

G4PDisplacedSolid::~G4PDisplacedSolid()
{;}

HepRef(G4PVSolid) G4PDisplacedSolid::GetConstituentMovedSolid()
{ return fPtrSolid; }

G4VSolid* G4PDisplacedSolid::MakeTransientObject() const
{ return 0; }

G4VSolid* G4PDisplacedSolid::MakeTransientDisplacedSolid
                                     (G4VSolid* movedSolid) const
{
  G4AffineTransform aTransform = fDirectTransform->MakeTransientObject();

  G4VSolid* theSolid = new G4DisplacedSolid
                              ( GetName(), movedSolid, aTransform );
  return theSolid;
}


