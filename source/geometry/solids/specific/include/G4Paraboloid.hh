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
// G4Paraboloid
//
// Class description:
//
//   A G4Paraboloid represents a solid with parabolic profile with possible
//   cuts along the Z axis.
//
//   Member Data:
//
//      dz              z half lenght
//      r1              radius at -dz
//      r2              radius at  dz
//      r2  > r1
//
//      Equation for the solid:
//        rho^2 <= k1 * z + k2;
//        -dz <= z <= dz
//        r1^2 = k1 * (-dz) + k2
//        r2^2 = k1 * ( dz) + k2

// 10.07.2007, L.Lindroos (CERN) - First implementation
// --------------------------------------------------------------------
#ifndef G4PARABOLOID_HH
#define G4PARABOLOID_HH

#include "G4GeomTypes.hh"

#if defined(G4GEOM_USE_USOLIDS)
#define G4GEOM_USE_UPARABOLOID 1
#endif

#if (defined(G4GEOM_USE_UPARABOLOID) && defined(G4GEOM_USE_SYS_USOLIDS))
  #define G4UParaboloid G4Paraboloid
  #include "G4UParaboloid.hh"
#else

#include <CLHEP/Units/PhysicalConstants.h>

#include "G4VSolid.hh"
#include "G4Polyhedron.hh"

class G4Paraboloid : public G4VSolid
{
  public:

    G4Paraboloid(const G4String& pName,
                       G4double  pDz,
                       G4double  pR1,
                       G4double  pR2);

    ~G4Paraboloid() override;

    // Access functions

    inline G4double GetZHalfLength() const;
    inline G4double GetRadiusMinusZ() const;
    inline G4double GetRadiusPlusZ() const;

    inline G4double GetCubicVolume() override;
    inline G4double GetSurfaceArea() override;
    inline G4double CalculateSurfaceArea() const;

    // Modifiers functions

    inline void SetZHalfLength(G4double dz);
    inline void SetRadiusMinusZ(G4double R1);
    inline void SetRadiusPlusZ(G4double R2);

    // Solid standard methods

    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const override;
    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                                 G4double& pmin, G4double& pmax) const override;
    EInside Inside(const G4ThreeVector& p) const override;
    G4ThreeVector SurfaceNormal( const G4ThreeVector& p) const override;
    G4double DistanceToIn(const G4ThreeVector& p,
                          const G4ThreeVector& v) const override;
    G4double DistanceToIn(const G4ThreeVector& p) const override;
    G4double DistanceToOut(const G4ThreeVector& p,
                           const G4ThreeVector& v,
                           const G4bool calcNorm = false,
                                 G4bool* validNorm = nullptr,
                                 G4ThreeVector* n = nullptr) const override;
    G4double DistanceToOut(const G4ThreeVector& p) const override;

    G4GeometryType GetEntityType() const override;

    G4VSolid* Clone() const override;

    std::ostream& StreamInfo(std::ostream& os) const override;

    G4ThreeVector GetPointOnSurface() const override;

    // Visualisation functions

    void DescribeYourselfTo(G4VGraphicsScene& scene) const override;
    G4Polyhedron* CreatePolyhedron() const override;
    G4Polyhedron* GetPolyhedron () const override;

    G4Paraboloid(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4Paraboloid(const G4Paraboloid& rhs);
    G4Paraboloid& operator=(const G4Paraboloid& rhs);
      // Copy constructor and assignment operator.

  protected:

    mutable G4bool fRebuildPolyhedron = false;
    mutable G4Polyhedron* fpPolyhedron = nullptr;

  private:

    // Making this mutable to allow GetPointOnSurface to have access to
    // area function.
    mutable G4double fSurfaceArea = 0.0;
    G4double fCubicVolume = 0.0;

    G4double dz, r1, r2;
    G4double k1, k2;
    // Defined to make some calculations easier to follow
};

#include "G4Paraboloid.icc"

#endif  // defined(G4GEOM_USE_UPARABOLOID) && defined(G4GEOM_USE_SYS_USOLIDS)

#endif // G4Paraboloid_HH
