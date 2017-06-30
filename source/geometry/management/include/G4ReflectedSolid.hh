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
// $Id: G4ReflectedSolid.hh 104317 2017-05-24 13:08:38Z gcosmo $
//
//
// class G4ReflectedSolid
//
// Class description:
//
// A Reflected solid is a solid that has been shifted from its original
// frame of reference to a new one.

// History:
//
// 23.07.01 V.Grichine: created
// --------------------------------------------------------------------
#ifndef G4ReflectedSolid_HH
#define G4ReflectedSolid_HH

#include "G4VSolid.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"

class G4ReflectedSolid : public G4VSolid
{
  public:  // with description

    G4ReflectedSolid( const G4String& pName,
                            G4VSolid* pSolid ,
                      const G4Transform3D& transform ) ;
      // For use in instantiating a transient instance.

    virtual ~G4ReflectedSolid();
      // Virtual destructor.

  public:  // without description 

    // Includes all the methods that a solid requires.

    EInside Inside( const G4ThreeVector& p ) const; 

    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const;

    G4bool CalculateExtent( const EAxis pAxis,
                            const G4VoxelLimits& pVoxelLimit,
                            const G4AffineTransform& pTransform,
                                  G4double& pMin, G4double& pMax) const;

    G4ThreeVector SurfaceNormal( const G4ThreeVector& p ) const;

    G4double DistanceToIn( const G4ThreeVector& p,
                           const G4ThreeVector& v ) const;

    G4double DistanceToIn( const G4ThreeVector& p) const;

    G4double DistanceToOut( const G4ThreeVector& p,
                            const G4ThreeVector& v,
                            const G4bool calcNorm=false,
                                  G4bool *validNorm=0,
                                  G4ThreeVector *n=0 ) const;

    G4double DistanceToOut( const G4ThreeVector& p ) const;

    void ComputeDimensions(       G4VPVParameterisation* p,
                            const G4int n,
                            const G4VPhysicalVolume* pRep );

    G4ThreeVector GetPointOnSurface() const;

    G4VSolid* Clone() const;

  public:  // with description 

    virtual G4GeometryType  GetEntityType() const;

    virtual const G4ReflectedSolid* GetReflectedSolidPtr() const;
    virtual       G4ReflectedSolid* GetReflectedSolidPtr();
      // If the Solid is a "G4ReflectedSolid",
      // return a self pointer else return 0.

    G4VSolid* GetConstituentMovedSolid() const;

    G4Transform3D GetTransform3D() const; 
    G4Transform3D GetDirectTransform3D() const; 
    void SetDirectTransform3D(G4Transform3D&);
      // Accessors methods.

    std::ostream& StreamInfo(std::ostream& os) const;

  public:  // without description

    G4ReflectedSolid(const G4ReflectedSolid& rhs);
    G4ReflectedSolid& operator=(const G4ReflectedSolid& rhs);
      // Copy constructor and assignment operator.

    void DescribeYourselfTo ( G4VGraphicsScene& scene ) const;
    G4Polyhedron* CreatePolyhedron () const;
    G4Polyhedron* GetPolyhedron    () const;
      // For creating graphical representations (i.e. for visualisation).

  protected:

    G4VSolid*          fPtrSolid;
    G4Transform3D*     fDirectTransform3D;

    mutable G4bool fRebuildPolyhedron;
    mutable G4Polyhedron* fpPolyhedron;  // Caches reflected G4Polyhedron.
};

#endif
