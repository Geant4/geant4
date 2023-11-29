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
// G4SubtractionSolid
//
// Class description:
//
// Class for description of subtraction of two solids: A - B.

// 14.10.98 V.Grichine: first implementation
// --------------------------------------------------------------------
#ifndef G4SUBTRACTIONSOLID_HH
#define G4SUBTRACTIONSOLID_HH

#include "G4BooleanSolid.hh"
#include "G4VSolid.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4AffineTransform.hh"


class G4SubtractionSolid : public G4BooleanSolid
{
  public:

    G4SubtractionSolid(  const G4String& pName,
                               G4VSolid* pSolidA ,
                               G4VSolid* pSolidB   ) ;

    G4SubtractionSolid(  const G4String& pName,
                               G4VSolid* pSolidA ,
                               G4VSolid* pSolidB ,
                               G4RotationMatrix* rotMatrix,
                         const G4ThreeVector& transVector   ) ;

    G4SubtractionSolid(  const G4String& pName,
                               G4VSolid* pSolidA ,
                               G4VSolid* pSolidB ,
                         const G4Transform3D& transform   ) ;

    ~G4SubtractionSolid() override ;

    G4GeometryType  GetEntityType() const override ;

    G4VSolid* Clone() const override;

    G4SubtractionSolid(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4SubtractionSolid(const G4SubtractionSolid& rhs);
    G4SubtractionSolid& operator=(const G4SubtractionSolid& rhs);
      // Copy constructor and assignment operator.

    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const override;

    G4bool CalculateExtent( const EAxis pAxis,
                            const G4VoxelLimits& pVoxelLimit,
                            const G4AffineTransform& pTransform,
                                  G4double& pMin, G4double& pMax) const override ;
       
    EInside Inside( const G4ThreeVector& p ) const override ;

    G4ThreeVector SurfaceNormal( const G4ThreeVector& p ) const override ;

    G4double DistanceToIn( const G4ThreeVector& p,
                           const G4ThreeVector& v  ) const override ;

    G4double DistanceToIn( const G4ThreeVector& p) const override ;

    G4double DistanceToOut( const G4ThreeVector& p,
                            const G4ThreeVector& v,
                            const G4bool calcNorm = false,
                                  G4bool* validNorm = nullptr,
                                  G4ThreeVector* n = nullptr ) const override ;

    G4double DistanceToOut( const G4ThreeVector& p ) const override ;


    void ComputeDimensions(       G4VPVParameterisation* p,
                            const G4int n,
                            const G4VPhysicalVolume* pRep ) override ;
                                   
    void DescribeYourselfTo ( G4VGraphicsScene& scene ) const override ;
    G4Polyhedron* CreatePolyhedron () const override ;

    G4double GetCubicVolume() final;
};

#endif
