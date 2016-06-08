// Persistent base class for solids created by boolean operations
// between other solids
//
// History:
// 10.11.99 Y.Morita, Initial Creation

#include "G4PBooleanSolid.hh"

G4PBooleanSolid::G4PBooleanSolid( const G4String& pName,
                                  HepRef(G4PVSolid) persSolidA,
                                  HepRef(G4PVSolid) persSolidB )
 : G4PVSolid(pName),
   createdDisplacedSolid(false)
{
   fPtrSolidA = persSolidA;
   fPtrSolidB = persSolidB;
}

G4PBooleanSolid::~G4PBooleanSolid()
{;}

///////////////////////////////////////////////////////////////
//
// If Solid is made up from a Boolean operation of two solids,
//   return the corresponding solid (for no=0 and 1)
// If the solid is not a "Boolean", return 0
const HepRef(G4PVSolid) G4PBooleanSolid::GetConstituentSolid(G4int no) const
{
  HepRef(G4PVSolid)  subSolid;
  if( no == 0 )
    subSolid = fPtrSolidA;
  else if( no == 1 )
    subSolid = fPtrSolidB;
  else
    G4Exception("G4PBooleanSolid::GetConstituentSolid()const invalid subsolid index");

  const HepRef(G4PVSolid) theSolid = subSolid;
  return theSolid;
}

  HepRef(G4PVSolid) G4PBooleanSolid::GetConstituentSolid(G4int no)
{
  HepRef(G4PVSolid)  subSolid;
  if( no == 0 )
    subSolid = fPtrSolidA;
  else if( no == 1 )
    subSolid = fPtrSolidB;
  else
    G4Exception("G4PBooleanSolid::GetConstituentSolid invalid subsolid index");

  return subSolid;
}

