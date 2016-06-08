// Persistent class for description of intersection of two CSG solids
//
// History: 
// 10.11.99 Y.Morita, Initial creation
//

#ifndef G4PUNIONSOLID_DDL
#define G4PUNIONSOLID_DDL

#include "G4PersistentSchema.hh"
#include "G4PBooleanSolid.hh"

class G4VSolid;

class G4PUnionSolid
 : public G4PBooleanSolid
{
public:
        G4PUnionSolid( const G4String& pName,
                       HepRef(G4PVSolid) persSolidA,
                       HepRef(G4PVSolid) persSolidB );
        virtual ~G4PUnionSolid() ;

        G4VSolid* MakeTransientObject() const;

        G4VSolid* MakeTransientBooleanSolid(
                               G4VSolid* aSolidA,
                               G4VSolid* aSolidB ) const;

        virtual G4GeometryType GetEntityType() const 
        { return G4String("G4UnionSolid"); }

protected:

private:

};

#endif

