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
// $Id: G4Box.hh 104316 2017-05-24 13:04:23Z gcosmo $
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// 
// G4Box
//
// Class description:
//
//   A Box is a cuboid of given half lengths dx,dy,dz. The Box is
//   centred on the origin with sides parallel to the x/y/z axes.

// History:
// 30.06.95 P.Kent: Converted from source code developed end 94
// 27.03.96 J.Allison: Added virtual functions DescribeYourselfTo() and
//                     SendWireframeTo(G4VGraphicsModel&)
// 22.07.96 J.Allison: Changed G4VGraphicsModel to G4VGraphicsScene
//                     and SendPolyhedronTo() to CreatePolyhedron()
// 27.03.98 J.Apostolakis: Inherit from G4CSGSolid (not G4VSolid)
// 18.11.99 J.Apostolakis, V.Grichine: kUndefined was added to ESide
// --------------------------------------------------------------------
#ifndef G4BOX_HH
#define G4BOX_HH

#if defined(G4GEOM_USE_USOLIDS)
#define G4GEOM_USE_UBOX 1
#endif

#if defined(G4GEOM_USE_UBOX)
  #define G4UBox G4Box
  #include "G4UBox.hh"
#else

#include "G4CSGSolid.hh"
#include "G4Polyhedron.hh"

class G4Box : public G4CSGSolid 
{
  public:  // with description

    G4Box(const G4String& pName, G4double pX, G4double pY, G4double pZ);
      // Construct a box with name, and half lengths pX,pY,pZ

    virtual ~G4Box();


    void ComputeDimensions(G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep);

    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const;

    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                                 G4double& pMin, G4double& pMax) const;

  // Accessors and modifiers

    inline G4double GetXHalfLength() const;
    inline G4double GetYHalfLength() const;
    inline G4double GetZHalfLength() const;

    void SetXHalfLength(G4double dx) ;
    void SetYHalfLength(G4double dy) ;
    void SetZHalfLength(G4double dz) ;

  // Methods for solid

    inline G4double GetCubicVolume();
    inline G4double GetSurfaceArea();

    EInside Inside(const G4ThreeVector& p) const;
    G4ThreeVector SurfaceNormal( const G4ThreeVector& p) const;
    G4double DistanceToIn(const G4ThreeVector& p, const G4ThreeVector& v) const;
    G4double DistanceToIn(const G4ThreeVector& p) const;
    G4double DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
                           const G4bool calcNorm=false,
                                 G4bool *validNorm=0, G4ThreeVector *n=0) const;
    G4double DistanceToOut(const G4ThreeVector& p) const;

    G4GeometryType GetEntityType() const;
    G4ThreeVector GetPointOnSurface() const; 

    G4VSolid* Clone() const;

    std::ostream& StreamInfo(std::ostream& os) const;

  // Utilities for visualization

    void          DescribeYourselfTo (G4VGraphicsScene& scene) const;
    G4VisExtent   GetExtent          () const;
    G4Polyhedron* CreatePolyhedron   () const;

  public:  // without description

    G4Box(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4Box(const G4Box& rhs);
    G4Box& operator=(const G4Box& rhs); 
      // Copy constructor and assignment operator.

  private:

    G4ThreeVector ApproxSurfaceNormal( const G4ThreeVector& p) const;
      // Algorithm for SurfaceNormal() following the original
      // specification for points not on the surface

  private:

    G4double fDx,fDy,fDz;
    G4double delta;  // Cached half Cartesian tolerance
};

#include "G4Box.icc"

#endif

#endif
