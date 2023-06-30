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
// G4UEllipticalCone
//
// Class description:
//
// Wrapper class for G4EllipticalCone to make use of VecGeom EllipticalCone.

// 13.09.19 Gabriele Cosmo, CERN
// --------------------------------------------------------------------
#ifndef G4UELLIPTICALCONE_HH
#define G4UELLIPTICALCONE_HH

#include "G4UAdapter.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include <VecGeom/volumes/UnplacedEllipticalCone.h>

#include "G4Polyhedron.hh"

class G4UEllipticalCone : public G4UAdapter<vecgeom::UnplacedEllipticalCone>
{
  using Shape_t = vecgeom::UnplacedEllipticalCone;
  using Base_t  = G4UAdapter<vecgeom::UnplacedEllipticalCone>;

  public:

    G4UEllipticalCone(const G4String& name, G4double pxSemiAxis,
                                            G4double pySemiAxis,
                                            G4double zMax,
                                            G4double pzTopCut);
   ~G4UEllipticalCone() override;

    G4VSolid* Clone() const override;

    G4double GetSemiAxisMin () const;
    G4double GetSemiAxisMax () const;
    G4double GetSemiAxisX () const;
    G4double GetSemiAxisY () const;
    G4double GetZMax() const;
    G4double GetZTopCut() const;
    void SetSemiAxis (G4double x, G4double y, G4double z);
    void SetZCut (G4double newzTopCut);

    inline G4GeometryType GetEntityType() const override;

    G4UEllipticalCone(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4UEllipticalCone( const G4UEllipticalCone& source );
    G4UEllipticalCone& operator=( const G4UEllipticalCone& source );
      // Copy constructor and assignment operator.

    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const override;
    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                           G4double& pmin, G4double& pmax) const override;
    G4Polyhedron* CreatePolyhedron() const override;
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline G4GeometryType G4UEllipticalCone::GetEntityType() const
{
  return "G4EllipticalCone";
}

#endif  // G4GEOM_USE_USOLIDS

#endif
