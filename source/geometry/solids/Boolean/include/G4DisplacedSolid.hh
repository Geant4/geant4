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
//
// $Id: G4DisplacedSolid.hh 104316 2017-05-24 13:04:23Z gcosmo $
//
//
// class G4DisplacedSolid
//
// Class description:
//
// A displaced solid is a solid that has been shifted from its original
// frame of reference to a new one. It is meant to be used only for
// simplifying the implementation of "Boolean solids". 

// History:
//
// 28.10.98 V.Grichine: created
// 22.11.00 V.Grichine: new set methods for matrix/vectors
//
// --------------------------------------------------------------------
#ifndef G4DisplacedSolid_HH
#define G4DisplacedSolid_HH

#include "G4VSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4AffineTransform.hh"

class G4DisplacedSolid : public G4VSolid
{
  public:  // with description

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

    virtual ~G4DisplacedSolid() ;

  public:  // without description 

    // It also has all the methods that a solid requires, eg.

    EInside Inside( const G4ThreeVector& p ) const ; 

    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const;

    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                                 G4double& pMin, G4double& pMax) const ;

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


    void ComputeDimensions(       G4VPVParameterisation* p,
                            const G4int n,
                            const G4VPhysicalVolume* pRep ) ;

    void CleanTransformations();

    G4ThreeVector GetPointOnSurface() const;

  public:  // with description 

    G4GeometryType  GetEntityType() const;
    G4VSolid* Clone() const;

    const G4DisplacedSolid* GetDisplacedSolidPtr() const;
          G4DisplacedSolid* GetDisplacedSolidPtr();
      // If the Solid is a "G4DisplacedSolid",
      // return a self pointer else return 0.

    G4VSolid*                GetConstituentMovedSolid() const;

    G4AffineTransform        GetTransform() const; 
    void       SetTransform(G4AffineTransform& ); 

    G4AffineTransform        GetDirectTransform() const; 
    void       SetDirectTransform(G4AffineTransform&); 
      // Access/Set methods.

    G4RotationMatrix         GetFrameRotation() const;
    void  SetFrameRotation(const G4RotationMatrix&);

    G4ThreeVector            GetFrameTranslation() const; 
    void  SetFrameTranslation(const G4ThreeVector&); 
      // Get/Set the rotation/translation, as applied to the frame of reference.

    G4RotationMatrix         GetObjectRotation() const;
    void  SetObjectRotation(const G4RotationMatrix&);

    G4ThreeVector            GetObjectTranslation() const; 
    void  SetObjectTranslation(const G4ThreeVector&); 
      // Get/Set the rotation/translation, as applied to the object.

    std::ostream& StreamInfo(std::ostream& os) const;

  public:  // without description

    G4DisplacedSolid(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4DisplacedSolid(const G4DisplacedSolid& rhs);
    G4DisplacedSolid& operator=(const G4DisplacedSolid& rhs);
      // Copy constructor and assignment operator.

    void DescribeYourselfTo ( G4VGraphicsScene& scene ) const ;
    G4Polyhedron* CreatePolyhedron () const ;
    G4Polyhedron* GetPolyhedron    () const ;
      // For creating graphical representations (ie for visualisation).

  protected:

    G4VSolid* fPtrSolid ;
    G4AffineTransform* fPtrTransform ;
    G4AffineTransform* fDirectTransform ;
    mutable G4bool fRebuildPolyhedron;
    mutable G4Polyhedron* fpPolyhedron;
} ;

#endif
