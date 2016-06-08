// Persistent class for description of subtraction of two CSG solids: A - B
//
// History: 
// 10.11.99 Y.Morita, Initial creation
//

#ifndef G4PSUBTRACTIONSOLID_DDL
#define G4PSUBTRACTIONSOLID_DDL

#include "G4PersistentSchema.hh"
#include "G4PBooleanSolid.hh"

class G4VSolid;

class G4PSubtractionSolid
 : public G4PBooleanSolid
{
public:
        G4PSubtractionSolid( const G4String& pName,
                             HepRef(G4PVSolid) persSolidA,
                             HepRef(G4PVSolid) persSolidB );
        virtual ~G4PSubtractionSolid() ;

        G4VSolid* MakeTransientObject() const;

        G4VSolid* MakeTransientBooleanSolid(
                               G4VSolid* aSolidA,
                               G4VSolid* aSolidB ) const;

        virtual G4GeometryType GetEntityType() const 
        { return G4String("G4SubtractionSolid"); }

protected:

private:

};

#endif

