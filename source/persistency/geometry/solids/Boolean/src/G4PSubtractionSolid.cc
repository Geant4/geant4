// Persistent class for description of subtraction of two CSG solids: A - B
//
// History: 
// 10.11.99 Y.Morita, Initial creation
//

#include "G4PSubtractionSolid.hh"
#include "G4SubtractionSolid.hh"

G4PSubtractionSolid::G4PSubtractionSolid( const G4String& pName,
                              HepRef(G4PVSolid) persSolidA,
                              HepRef(G4PVSolid) persSolidB )
 : G4PBooleanSolid(pName, persSolidA, persSolidB)
{;}

G4PSubtractionSolid::~G4PSubtractionSolid()
{;}

G4VSolid* G4PSubtractionSolid::MakeTransientObject() const
{ return 0; }

G4VSolid* G4PSubtractionSolid::MakeTransientBooleanSolid(
                              G4VSolid* aSolidA,
                              G4VSolid* aSolidB ) const
{
  return new G4SubtractionSolid( GetName(), aSolidA, aSolidB );
}


