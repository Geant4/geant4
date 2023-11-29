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
// G4ReflectedSolid
//
// Class description:
//
// A Reflected solid is a solid that has been shifted from its original
// frame of reference to a new one.

// 23.07.01, V.Grichine - created
// --------------------------------------------------------------------
#ifndef G4ReflectedSolid_HH
#define G4ReflectedSolid_HH

#include "G4VSolid.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"

class G4ReflectedSolid : public G4VSolid
{
  public:

    G4ReflectedSolid( const G4String& pName,
                            G4VSolid* pSolid ,
                      const G4Transform3D& transform ) ;
      // For use in instantiating a transient instance.

    ~G4ReflectedSolid() override;
      // Virtual destructor.

    // Includes all the methods that a solid requires.

    EInside Inside( const G4ThreeVector& p ) const override; 

    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const override;

    G4bool CalculateExtent( const EAxis pAxis,
                            const G4VoxelLimits& pVoxelLimit,
                            const G4AffineTransform& pTransform,
                                  G4double& pMin, G4double& pMax) const override;

    G4ThreeVector SurfaceNormal( const G4ThreeVector& p ) const override;

    G4double DistanceToIn( const G4ThreeVector& p,
                           const G4ThreeVector& v ) const override;

    G4double DistanceToIn( const G4ThreeVector& p) const override;

    G4double DistanceToOut( const G4ThreeVector& p,
                            const G4ThreeVector& v,
                            const G4bool calcNorm = false,
                                  G4bool* validNorm = nullptr,
                                  G4ThreeVector* n = nullptr ) const override;

    G4double DistanceToOut( const G4ThreeVector& p ) const override;

    void ComputeDimensions(       G4VPVParameterisation* p,
                            const G4int n,
                            const G4VPhysicalVolume* pRep ) override;

    G4double GetCubicVolume() override;
    G4double GetSurfaceArea() override;

    G4ThreeVector GetPointOnSurface() const override;

    G4VSolid* Clone() const override;

    G4GeometryType  GetEntityType() const override;

    virtual const G4ReflectedSolid* GetReflectedSolidPtr() const;
    virtual       G4ReflectedSolid* GetReflectedSolidPtr();
      // If the Solid is a "G4ReflectedSolid",
      // return a self pointer else return nullptr.

    G4VSolid* GetConstituentMovedSolid() const;

    G4Transform3D GetTransform3D() const; 
    G4Transform3D GetDirectTransform3D() const; 
    void SetDirectTransform3D(G4Transform3D&);
      // Accessors methods.

    std::ostream& StreamInfo(std::ostream& os) const override;

    G4ReflectedSolid(const G4ReflectedSolid& rhs);
    G4ReflectedSolid& operator=(const G4ReflectedSolid& rhs);
      // Copy constructor and assignment operator.

    void DescribeYourselfTo ( G4VGraphicsScene& scene ) const override;
    G4Polyhedron* CreatePolyhedron () const override;
    G4Polyhedron* GetPolyhedron    () const override;
      // For creating graphical representations (i.e. for visualisation).

  protected:

    G4VSolid*          fPtrSolid = nullptr;
    G4Transform3D*     fDirectTransform3D = nullptr;

    mutable G4bool fRebuildPolyhedron = false;
    mutable G4Polyhedron* fpPolyhedron = nullptr;
      // Caches reflected G4Polyhedron.
};

#endif
