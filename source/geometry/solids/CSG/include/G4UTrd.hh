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
// G4UTrd
//
// Class description:
//
// Wrapper class for G4Trd to make use of VecGeom Trd.

// 13.09.13 G.Cosmo, CERN
// --------------------------------------------------------------------
#ifndef G4UTRD_HH
#define G4UTRD_HH

#include "G4UAdapter.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include <VecGeom/volumes/UnplacedTrd.h>

#include "G4Polyhedron.hh"

class G4UTrd : public G4UAdapter<vecgeom::GenericUnplacedTrd> 
{
  using Shape_t = vecgeom::GenericUnplacedTrd;
  using Base_t  = G4UAdapter<vecgeom::GenericUnplacedTrd>;

  public:

    G4UTrd(const G4String& pName,
                 G4double pdx1, G4double pdx2,
                 G4double pdy1, G4double pdy2,
                 G4double pdz);
      // Constructs a trapezoid with name, and half lengths

   ~G4UTrd() override;

    void ComputeDimensions(      G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep) override;

    G4VSolid* Clone() const override;

    G4double GetXHalfLength1() const;
    G4double GetXHalfLength2() const;
    G4double GetYHalfLength1() const;
    G4double GetYHalfLength2() const;
    G4double GetZHalfLength()  const;

    void SetXHalfLength1(G4double val);
    void SetXHalfLength2(G4double val);
    void SetYHalfLength1(G4double val);
    void SetYHalfLength2(G4double val);
    void SetZHalfLength(G4double val);

    void SetAllParameters(G4double pdx1, G4double pdx2,
                          G4double pdy1, G4double pdy2, G4double pdz);

    inline G4GeometryType GetEntityType() const override;

    inline G4bool IsFaceted() const override;

    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const override;

    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                           G4double& pMin, G4double& pMax) const override;

    G4Polyhedron* CreatePolyhedron() const override;

    G4UTrd(const G4UTrd& rhs);
    G4UTrd& operator=(const G4UTrd& rhs); 
      // Copy constructor and assignment operator.
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline G4GeometryType G4UTrd::GetEntityType() const
{
  return "G4Trd";
}

inline G4bool G4UTrd::IsFaceted() const
{
  return true;
}

#endif  // G4GEOM_USE_USOLIDS

#endif
