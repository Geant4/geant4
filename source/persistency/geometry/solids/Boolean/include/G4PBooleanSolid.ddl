// Persistent base class for solids created by boolean operations
// between other solids
//
// History:
// 10.11.99 Y.Morita, Initial Creation

#ifndef G4PBOOLEANSOLID_DDL
#define G4PBOOLEANSOLID_DDL

#include "G4PVSolid.hh"
#include "G4PersistentSchema.hh"

class G4PBooleanSolid
 : public G4PVSolid
{
public:
          G4PBooleanSolid( const G4String& pName,
                           HepRef(G4PVSolid) persSolidA,
                           HepRef(G4PVSolid) persSolidB );
          virtual ~G4PBooleanSolid();

          virtual G4VSolid* MakeTransientBooleanSolid(
                               G4VSolid* aSolidA,
                               G4VSolid* aSolidB ) const = 0;

public:  // With Description
    // If Solid is made up from a Boolean operation of two solids,
    //   return the corresponding solid (for no=0 and 1)
    // If the solid is not a "Boolean", return 0
    virtual const HepRef(G4PVSolid) GetConstituentSolid(G4int no) const;
    virtual       HepRef(G4PVSolid) GetConstituentSolid(G4int no);

protected:
          d_Ref<G4PVSolid> fPtrSolidA;
          d_Ref<G4PVSolid> fPtrSolidB;

private:
          G4bool  createdDisplacedSolid;
};

#endif

