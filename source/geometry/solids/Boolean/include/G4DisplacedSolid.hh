//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// G4DisplacedSolid
//
// Class description:
//
// A displaced solid is a solid that has been shifted from its original
// frame of reference to a new one. It is meant to be used only for
// simplifying the implementation of "Boolean solids". 

// 28.10.98 V.Grichine - created
// --------------------------------------------------------------------
#ifndef G4DISPLACEDSOLID_HH
#define G4DISPLACEDSOLID_HH

#include "G4VSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4AffineTransform.hh"

class G4DisplacedSolid : public G4VSolid
{
  public:

    G4DisplacedSolid( const G4String& pName,
                            G4VSolid* pSolid ,
                            G4RotationMatrix* rotMatrix,
                      const G4ThreeVector& transVector  ) ;

    G4DisplacedSolid( const G4String& pName,
                            G4VSolid* pSolid ,
                      const G4Transform3D& transform  ) ;

    G4DisplacedSolid( const G4String& pName,
                            G4VSolid* pSolid ,
                      const G4AffineTransform directTransform );
      // For use in instantiating a transient instance from a persistent one.

    ~G4DisplacedSolid() override ;

    EInside Inside( const G4ThreeVector& p ) const override ; 

    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const override;

    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                                 G4double& pMin, G4double& pMax) const override ;

    G4ThreeVector SurfaceNormal( const G4ThreeVector& p ) const override ;

    G4double DistanceToIn( const G4ThreeVector& p,
                           const G4ThreeVector& v  ) const override ;

    G4double DistanceToIn( const G4ThreeVector& p) const override ;

    G4double DistanceToOut( const G4ThreeVector& p,
                            const G4ThreeVector& v,
                            const G4bool calcNorm=false,
                                  G4bool *validNorm=nullptr,
                                  G4ThreeVector *n=nullptr ) const override ;

    G4double DistanceToOut( const G4ThreeVector& p ) const override ;


    void ComputeDimensions(       G4VPVParameterisation* p,
                            const G4int n,
                            const G4VPhysicalVolume* pRep ) override ;

    void CleanTransformations();

    G4double GetCubicVolume() override;
    G4double GetSurfaceArea() override;

    G4ThreeVector GetPointOnSurface() const override;

    G4GeometryType  GetEntityType() const override;
    G4VSolid* Clone() const override;

    const G4DisplacedSolid* GetDisplacedSolidPtr() const override;
          G4DisplacedSolid* GetDisplacedSolidPtr() override;
      // If the Solid is a "G4DisplacedSolid",
      // return a self pointer else return 0.

    G4VSolid* GetConstituentMovedSolid() const;

    G4AffineTransform GetTransform() const; 
    void       SetTransform(G4AffineTransform& ); 

    G4AffineTransform GetDirectTransform() const; 
    void       SetDirectTransform(G4AffineTransform&); 
      // Access/Set methods.

    G4RotationMatrix GetFrameRotation() const;
    void SetFrameRotation(const G4RotationMatrix&);

    G4ThreeVector GetFrameTranslation() const; 
    void SetFrameTranslation(const G4ThreeVector&); 
      // Get/Set the rotation/translation, as applied to the frame of reference.

    G4RotationMatrix GetObjectRotation() const;
    void SetObjectRotation(const G4RotationMatrix&);

    G4ThreeVector GetObjectTranslation() const; 
    void SetObjectTranslation(const G4ThreeVector&); 
      // Get/Set the rotation/translation, as applied to the object.

    std::ostream& StreamInfo(std::ostream& os) const override;

    G4DisplacedSolid(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4DisplacedSolid(const G4DisplacedSolid& rhs);
    G4DisplacedSolid& operator=(const G4DisplacedSolid& rhs);
      // Copy constructor and assignment operator.

    void DescribeYourselfTo ( G4VGraphicsScene& scene ) const override ;
    G4Polyhedron* CreatePolyhedron () const override ;
    G4Polyhedron* GetPolyhedron () const override ;
      // For creating graphical representations (ie for visualisation).

  protected:

    G4VSolid* fPtrSolid = nullptr;
    G4AffineTransform* fPtrTransform = nullptr;
    G4AffineTransform* fDirectTransform = nullptr;
    mutable G4bool fRebuildPolyhedron = false;
    mutable G4Polyhedron* fpPolyhedron = nullptr;
};

#endif
