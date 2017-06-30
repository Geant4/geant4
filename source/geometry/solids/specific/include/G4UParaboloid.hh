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
//
// G4UParaboloid
//
// Class description:
//
//   Wrapper class for UParaboloid to make use of it from USolids module.

// History:
// 19.08.15 Guilherme Lima, FNAL
// --------------------------------------------------------------------
#ifndef G4UPARABOLOID_HH
#define G4UPARABOLOID_HH

#include "G4UAdapter.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include <volumes/UnplacedParaboloid.h>

#include "G4Polyhedron.hh"

class G4UParaboloid : public G4UAdapter<vecgeom::UnplacedParaboloid>
{
  using Shape_t = vecgeom::UnplacedParaboloid;
  using Base_t  = G4UAdapter<vecgeom::UnplacedParaboloid>;

  public:  // with description

    G4UParaboloid(const G4String& name, G4double dz,
                                        G4double rlo,
                                        G4double rhi);
   ~G4UParaboloid();

    G4VSolid* Clone() const;

    G4double GetZHalfLength() const;
    G4double GetRadiusMinusZ() const;
    G4double GetRadiusPlusZ() const;

    void SetZHalfLength(G4double dz);
    void SetRadiusMinusZ(G4double r1);
    void SetRadiusPlusZ(G4double r2);

    inline G4GeometryType GetEntityType() const;

  public:  // without description

    G4UParaboloid(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4UParaboloid( const G4UParaboloid &source );
    G4UParaboloid &operator=( const G4UParaboloid &source );
      // Copy constructor and assignment operator.

    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const;
    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                           G4double& pmin, G4double& pmax) const;
    G4Polyhedron* CreatePolyhedron() const;
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline G4GeometryType G4UParaboloid::GetEntityType() const
{
  return "G4Paraboloid";
}

#endif  // G4GEOM_USE_USOLIDS

#endif
