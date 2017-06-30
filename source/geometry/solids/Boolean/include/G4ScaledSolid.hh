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
// $Id:$
//
//
// class G4ScaledSolid
//
// Class description:
//
// A scaled solid is a solid that has been scaled in dimensions
// in X, Y or Z, from its original description.

// History:
//
// 27.10.15 G.Cosmo: created
//
// --------------------------------------------------------------------
#ifndef G4ScaledSolid_HH
#define G4ScaledSolid_HH

#include "G4VSolid.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4AffineTransform.hh"

class G4ScaleTransform;

class G4ScaledSolid : public G4VSolid
{
  public:  // with description

    G4ScaledSolid( const G4String& pName,
                         G4VSolid* pSolid ,
                   const G4Scale3D& pScale  );

    virtual ~G4ScaledSolid();

    EInside Inside( const G4ThreeVector& p ) const;

    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const;

    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                                 G4double& pMin, G4double& pMax) const;

    G4ThreeVector SurfaceNormal( const G4ThreeVector& p ) const;

    G4double DistanceToIn( const G4ThreeVector& p,
                           const G4ThreeVector& v  ) const;

    G4double DistanceToIn( const G4ThreeVector& p) const;

    G4double DistanceToOut( const G4ThreeVector& p,
                            const G4ThreeVector& v,
                            const G4bool calcNorm=false,
                                  G4bool *validNorm=0,
                                  G4ThreeVector *n=0      ) const;

    G4double DistanceToOut( const G4ThreeVector& p ) const;

    void ComputeDimensions(       G4VPVParameterisation* p,
                            const G4int n,
                            const G4VPhysicalVolume* pRep );

    void CleanTransformations();

    G4ThreeVector GetPointOnSurface() const;

    G4Scale3D GetScaleTransform() const; 
    void SetScaleTransform(const G4Scale3D& scale); 

    G4VSolid* GetUnscaledSolid() const;

    G4GeometryType  GetEntityType() const;
    G4VSolid* Clone() const;

    std::ostream& StreamInfo(std::ostream& os) const;

  public:  // without description

    G4ScaledSolid(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4ScaledSolid(const G4ScaledSolid& rhs);
    G4ScaledSolid& operator=(const G4ScaledSolid& rhs);
      // Copy constructor and assignment operator.

    void DescribeYourselfTo ( G4VGraphicsScene& scene ) const;
    G4Polyhedron* CreatePolyhedron () const;
    G4Polyhedron* GetPolyhedron    () const;
      // For creating graphical representations (i.e. for visualisation).

  private:

    G4VSolid* fPtrSolid;
    G4ScaleTransform* fScale;
    mutable G4bool fRebuildPolyhedron;
    mutable G4Polyhedron* fpPolyhedron;
} ;

#endif
