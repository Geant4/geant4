// Class for description of intersection of two CSG solids
//
// History: 
//
// 12.09.98 V.Grichine
//

#ifndef G4UNIONSOLID_HH
#define G4UNIONSOLID_HH

#include "G4BooleanSolid.hh"
#include "G4VSolid.hh"
// #include "G4PlacedSolid.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4AffineTransform.hh"

class G4UnionSolid: public G4BooleanSolid
{
public:
                  G4UnionSolid( const G4String& pName,
                                       G4VSolid* pSolidA ,
                                       G4VSolid* pSolidB   ) ;

                  G4UnionSolid( const G4String& pName,
                                      G4VSolid* pSolidA ,
                                      G4VSolid* pSolidB ,
                                      G4RotationMatrix* rotMatrix,
                                const G4ThreeVector& transVector    ) ;

                  G4UnionSolid( const G4String& pName,
                                      G4VSolid* pSolidA ,
                                      G4VSolid* pSolidB ,
                                const G4Transform3D& transform   ) ;

                 virtual ~G4UnionSolid() ;


     G4bool CalculateExtent(const EAxis pAxis,
				   const G4VoxelLimits& pVoxelLimit,
				   const G4AffineTransform& pTransform,
				   G4double& pMin, G4double& pMax) const ;
       
    EInside Inside( const G4ThreeVector& p ) const ;

    G4ThreeVector SurfaceNormal( const G4ThreeVector& p ) const ;

    G4double DistanceToIn( const G4ThreeVector& p,
                           const G4ThreeVector& v  ) const ;

    G4double DistanceToIn( const G4ThreeVector& p) const ;

    G4double DistanceToOut( const G4ThreeVector& p,
			    const G4ThreeVector& v,
			    const G4bool calcNorm=false,
			    G4bool *validNorm=0,
			    G4ThreeVector *n=0      ) const ;

    G4double DistanceToOut( const G4ThreeVector& p ) const ;


    void ComputeDimensions( G4VPVParameterisation* p,
	                    const G4int n,
                            const G4VPhysicalVolume* pRep ) ;
                                   
    virtual G4GeometryType  GetEntityType() const 
    { return G4String("G4UnionSolid"); }

    void DescribeYourselfTo ( G4VGraphicsScene& scene ) const ;
    G4Polyhedron* CreatePolyhedron () const ;
    G4NURBS*      CreateNURBS      () const ;

protected:

private:

} ;

#endif

