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
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// G4UCons
//
// Class description:
//
//   Wrapper class for UCons to make use of UCons from USolids module.

// History:
// 30.10.13 G.Cosmo, CERN/PH
// --------------------------------------------------------------------

#ifndef G4UCONS_HH
#define G4UCONS_HH

#include "G4UAdapter.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include <volumes/UnplacedCone.h>

#include "G4Polyhedron.hh"

class G4UCons : public G4UAdapter<vecgeom::UnplacedCone>
{
  using Shape_t = vecgeom::UnplacedCone;
  using Base_t = G4UAdapter<vecgeom::UnplacedCone>;

  public:  // with description

    G4UCons(const G4String& pName,
                  G4double pRmin1, G4double pRmax1,
                  G4double pRmin2, G4double pRmax2,
                  G4double pDz,
                  G4double pSPhi, G4double pDPhi);
      // Constructs a cone with the given name and dimensions

   ~G4UCons();

    void ComputeDimensions(       G4VPVParameterisation* p,
                            const G4int n,
                            const G4VPhysicalVolume* pRep );

    G4VSolid* Clone() const;

    G4double GetInnerRadiusMinusZ() const;
    G4double GetOuterRadiusMinusZ() const;
    G4double GetInnerRadiusPlusZ()  const;
    G4double GetOuterRadiusPlusZ()  const;
    G4double GetZHalfLength()       const;
    G4double GetStartPhiAngle()     const;
    G4double GetDeltaPhiAngle()     const;
    G4double GetSinStartPhi()       const;
    G4double GetCosStartPhi()       const;
    G4double GetSinEndPhi()         const;
    G4double GetCosEndPhi()         const;
  
    void SetInnerRadiusMinusZ (G4double Rmin1 );
    void SetOuterRadiusMinusZ (G4double Rmax1 );
    void SetInnerRadiusPlusZ  (G4double Rmin2 );
    void SetOuterRadiusPlusZ  (G4double Rmax2 );
    void SetZHalfLength       (G4double newDz );
    void SetStartPhiAngle     (G4double newSPhi, G4bool trig=true);
    void SetDeltaPhiAngle     (G4double newDPhi);

    inline G4GeometryType GetEntityType() const;

    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const;

    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                                 G4double& pMin, G4double& pMax) const;

    G4Polyhedron* CreatePolyhedron() const;

  public:  // without description
       
    G4UCons(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4UCons(const G4UCons& rhs);
    G4UCons& operator=(const G4UCons& rhs); 
      // Copy constructor and assignment operator.
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline G4GeometryType G4UCons::GetEntityType() const
{
  return "G4Cons";
}

#endif  // G4GEOM_USE_USOLIDS

#endif
