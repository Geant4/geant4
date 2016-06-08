// Persistent class for description of intersection of two CSG solids
//
// History: 
// 10.11.99 Y.Morita, Initial creation
//

#include "G4PUnionSolid.hh"
#include "G4UnionSolid.hh"

G4PUnionSolid::G4PUnionSolid( const G4String& pName,
                              HepRef(G4PVSolid) persSolidA,
                              HepRef(G4PVSolid) persSolidB )
 : G4PBooleanSolid(pName, persSolidA, persSolidB)
{;}

G4PUnionSolid::~G4PUnionSolid()
{;}

G4VSolid* G4PUnionSolid::MakeTransientObject() const
{ return 0; }

G4VSolid* G4PUnionSolid::MakeTransientBooleanSolid(
                              G4VSolid* aSolidA,
                              G4VSolid* aSolidB ) const
{
  return new G4UnionSolid( GetName(), aSolidA, aSolidB );
}

