// Persistent class for description of intersection of two CSG solids
//
// History: 
// 10.11.99 Y.Morita, Initial creation

#include "G4PIntersectionSolid.hh"
#include "G4IntersectionSolid.hh"


G4PIntersectionSolid::G4PIntersectionSolid
                            ( const G4String& pName,
                              HepRef(G4PVSolid) persSolidA,
                              HepRef(G4PVSolid) persSolidB )
 : G4PBooleanSolid(pName, persSolidA, persSolidB)
{;}

G4PIntersectionSolid::~G4PIntersectionSolid()
{;}

G4VSolid* G4PIntersectionSolid::MakeTransientObject() const
{ return 0; }

G4VSolid* G4PIntersectionSolid::MakeTransientBooleanSolid(
                              G4VSolid* aSolidA,
                              G4VSolid* aSolidB ) const
{
  return new G4IntersectionSolid( GetName(), aSolidA, aSolidB );
}


