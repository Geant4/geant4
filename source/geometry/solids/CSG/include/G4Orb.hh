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
// G4Orb
//
// Class description:
//
//   A G4Orb is a simple case of G4Sphere. It has only:
//   fRmax  outer radius

// 20.08.03 V.Grichine - created
// 08.08.17 E.Tcherniaev - revised
// --------------------------------------------------------------------
#ifndef G4ORB_HH
#define G4ORB_HH

#include "G4GeomTypes.hh"

#if defined(G4GEOM_USE_USOLIDS)
#define G4GEOM_USE_UORB 1
#endif

#if defined(G4GEOM_USE_UORB)
  #define G4UOrb G4Orb
  #include "G4UOrb.hh"
#else

#include <CLHEP/Units/PhysicalConstants.h>

#include "G4CSGSolid.hh"
#include "G4Polyhedron.hh"

class G4Orb : public G4CSGSolid
{
  public:  // with description

    G4Orb(const G4String& pName, G4double pRmax);

    ~G4Orb();

  // Accessors and modifiers

    inline G4double GetRadius() const;
    inline G4double GetRadialTolerance() const;

    inline void SetRadius(G4double newRmax);

  // Methods for solid

    inline G4double GetCubicVolume();
    inline G4double GetSurfaceArea();

    void ComputeDimensions(      G4VPVParameterisation* p,
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

    G4ThreeVector GetPointOnSurface() const;

    G4VSolid* Clone() const;

    std::ostream& StreamInfo(std::ostream& os) const;

    // Visualisation functions

    void          DescribeYourselfTo (G4VGraphicsScene& scene) const;
    G4VisExtent   GetExtent          () const;
    G4Polyhedron* CreatePolyhedron   () const;

  public:  // without description

    G4Orb(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects

    G4Orb(const G4Orb& rhs);
    G4Orb& operator=(const G4Orb& rhs);
      // Copy constructor and assignment operator

  protected:

    void Initialize();

  private:

    G4double fRmax = 0.0;
    G4double halfRmaxTol = 0.0;
    G4double sqrRmaxPlusTol = 0.0, sqrRmaxMinusTol = 0.0;
};

#include "G4Orb.icc"

#endif

#endif
