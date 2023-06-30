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
// G4UEllipticalTube
//
// Class description:
//
// Wrapper class for G4EllipticalTube to make use of VecGeom EllipticalTube.

// 13.09.19 Gabriele Cosmo, CERN
// --------------------------------------------------------------------
#ifndef G4UELLIPTICALTUBE_HH
#define G4UELLIPTICALTUBE_HH

#include "G4UAdapter.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include <VecGeom/volumes/UnplacedEllipticalTube.h>

#include "G4Polyhedron.hh"

class G4UEllipticalTube : public G4UAdapter<vecgeom::UnplacedEllipticalTube>
{
  using Shape_t = vecgeom::UnplacedEllipticalTube;
  using Base_t  = G4UAdapter<vecgeom::UnplacedEllipticalTube>;

  public:

    G4UEllipticalTube(const G4String& name, G4double dx,
                                            G4double dy,
                                            G4double dz);
   ~G4UEllipticalTube() override;

    G4VSolid* Clone() const override;

    G4double GetDx() const;
    G4double GetDy() const;
    G4double GetDz() const;

    void SetDx(G4double dx);
    void SetDy(G4double dy);
    void SetDz(G4double dz);

    inline G4GeometryType GetEntityType() const override;

    G4UEllipticalTube(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4UEllipticalTube( const G4UEllipticalTube& source );
    G4UEllipticalTube &operator=( const G4UEllipticalTube& source );
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

inline G4GeometryType G4UEllipticalTube::GetEntityType() const
{
  return "G4EllipticalTube";
}

#endif  // G4GEOM_USE_USOLIDS

#endif
