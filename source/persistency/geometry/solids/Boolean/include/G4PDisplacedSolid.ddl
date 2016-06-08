// Persistent class describing solid placements for boolean operations
//
// History:
// 10.11.99 Y.Morita, Initial creation

#ifndef G4PDisplacedSolid_DDL
#define G4PDisplacedSolid_DDL

#include "G4PersistentSchema.hh"
#include "G4PVSolid.hh"
#include "G4PAffineTransform.hh"

class G4VSolid;

class G4PDisplacedSolid
 : public G4PVSolid
{
public:
        G4PDisplacedSolid ( HepRef(G4PVSolid) persCostituentSolid,
                            HepRef(G4PAffineTransform) pDirectTransform );
        virtual ~G4PDisplacedSolid();

        G4VSolid* MakeTransientObject() const;

        G4VSolid* MakeTransientDisplacedSolid(G4VSolid* aSolid) const;

        virtual G4GeometryType  GetEntityType() const
        { return G4String("G4DisplacedSolid"); }

        HepRef(G4PVSolid) GetConstituentMovedSolid();

protected:
        d_Ref<G4PVSolid> fPtrSolid;
//      fPtrTransform can be created by fDirectTransform.Inverse()
        d_Ref<G4PAffineTransform> fDirectTransform;

private:

};

#endif
