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
// G4Ellipsoid
//
// Class description:
//
//   A G4Ellipsoid is an ellipsoidal solid, optionally cut at a given z
//
//   Member Data:
//
//     xSemiAxis   semi-axis, X
//     ySemiAxis   semi-axis, Y
//     zSemiAxis   semi-axis, Z
//     zBottomCut  lower cut in Z (solid lies above this plane)
//     zTopCut     upper cut in Z (solid lies below this plane)

// 10.11.1999  G.Horton-Smith (Caltech, USA) - First implementation
// 10.02.2005  G.Guerrieri (INFN Genova, Italy) - Revision
// 15.12.2019  E.Tcherniaev - Complete revision
// --------------------------------------------------------------------
#ifndef G4ELLIPSOID_HH
#define G4ELLIPSOID_HH

#include "G4GeomTypes.hh"

#if defined(G4GEOM_USE_USOLIDS)
#define G4GEOM_USE_UELLIPSOID 1
#endif

#if (defined(G4GEOM_USE_UELLIPSOID) && defined(G4GEOM_USE_SYS_USOLIDS))
  #define G4UEllipsoid G4Ellipsoid
  #include "G4UEllipsoid.hh"
#else

#include <CLHEP/Units/PhysicalConstants.h>

#include "G4VSolid.hh"
#include "G4Polyhedron.hh"

class G4Ellipsoid : public G4VSolid
{
  public:

    // Constructor 
    G4Ellipsoid(const G4String& name,
                      G4double  xSemiAxis,
                      G4double  ySemiAxis,
                      G4double  zSemiAxis,
                      G4double  zBottomCut = 0.,
                      G4double  zTopCut = 0.);

    // Destructor
    virtual ~G4Ellipsoid();

    // Accessors
    inline G4double GetDx() const;
    inline G4double GetDy() const;
    inline G4double GetDz() const;
    inline G4double GetSemiAxisMax (G4int i) const;
    inline G4double GetZBottomCut() const;
    inline G4double GetZTopCut() const;

    // Modifiers
    inline void SetSemiAxis (G4double x, G4double y, G4double z);
    inline void SetZCuts (G4double newzBottomCut, G4double newzTopCut);

    // Standard methods
    void ComputeDimensions(G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep);

    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const;
    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                                 G4double& pmin, G4double& pmax) const;

    EInside Inside(const G4ThreeVector& p) const;
    G4ThreeVector SurfaceNormal( const G4ThreeVector& p) const;
    G4double DistanceToIn(const G4ThreeVector& p,
                          const G4ThreeVector& v) const;
    G4double DistanceToIn(const G4ThreeVector& p) const;
    G4double DistanceToOut(const G4ThreeVector& p,
                           const G4ThreeVector& v,
                           const G4bool calcNorm = false,
                                 G4bool* validNorm = nullptr,
                                 G4ThreeVector* n = nullptr) const;
    G4double DistanceToOut(const G4ThreeVector& p) const;

    G4GeometryType GetEntityType() const;

    G4VSolid* Clone() const;

    std::ostream& StreamInfo(std::ostream& os) const;

    G4double GetCubicVolume();
    G4double GetSurfaceArea();

    G4ThreeVector GetPointOnSurface() const;

    // Visualisation methods
    void DescribeYourselfTo(G4VGraphicsScene& scene) const;
    G4VisExtent GetExtent() const;
    G4Polyhedron* CreatePolyhedron() const;
    G4Polyhedron* GetPolyhedron () const;

  public:

    // Fake default constructor for usage restricted to direct object
    // persistency for clients requiring preallocation of memory for
    // persistifiable objects
    G4Ellipsoid(__void__&);

    // Copy constructorassignment operator
    G4Ellipsoid(const G4Ellipsoid& rhs);

    // Assignment
    G4Ellipsoid& operator=(const G4Ellipsoid& rhs);

  private:

    // Check parameters and set cached values
    void CheckParameters();

    // Return normal to surface closest to p
    G4ThreeVector ApproxSurfaceNormal(const G4ThreeVector& p) const;

    // Calculate area of lateral surface
    G4double LateralSurfaceArea() const;

  private:

    // Ellipsoid parameters
    G4double fDx;         // X semi-axis
    G4double fDy;         // Y semi-axis
    G4double fDz;         // Z semi-axis
    G4double fZBottomCut; // Bottom cut in Z
    G4double fZTopCut;    // Top cut in Z

    // Precalculated cached values
    G4double halfTolerance; // Surface tolerance
    G4double fXmax;         // X extent
    G4double fYmax;         // Y extent
    G4double fRsph;         // Radius of bounding sphere
    G4double fR;            // Radius of sphere after scaling
    G4double fSx;           // X scale factor
    G4double fSy;           // Y scale factor
    G4double fSz;           // Z scale factor
    G4double fZMidCut;      // Middle position between cuts after scaling
    G4double fZDimCut;      // Half distance between cut after scaling
    G4double fQ1; // 1st coefficient in approximation of dist = Q1*(x^2+y^2+z^2) - Q2
    G4double fQ2; // 2nd coefficient in approximation of dist = Q1*(x^2+y^2+z^2) - Q2

    G4double fCubicVolume = 0.0; // Volume
    G4double fSurfaceArea = 0.0; // Surface area
    mutable G4double fLateralArea = 0.0; // Lateral surface area
    mutable G4bool fRebuildPolyhedron = false;
    mutable G4Polyhedron* fpPolyhedron = nullptr;
};

#include "G4Ellipsoid.icc"

#endif  // defined(G4GEOM_USE_UELLIPSOID) && defined(G4GEOM_USE_SYS_USOLIDS)

#endif // G4ELLIPSOID_HH
