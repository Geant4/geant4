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
// G4UTet
//
// Class description:
//
// Wrapper class for G4Tet to make use of VecGeom Tet.

// 1.11.13 G.Cosmo, CERN
// --------------------------------------------------------------------
#ifndef G4UTET_HH
#define G4UTET_HH

#include "G4UAdapter.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include <VecGeom/volumes/UnplacedTet.h>

#include "G4Polyhedron.hh"

class G4UTet : public G4UAdapter<vecgeom::UnplacedTet>
{

  using Shape_t = vecgeom::UnplacedTet;
  using Base_t = G4UAdapter<vecgeom::UnplacedTet>;

  public:  // with description

    G4UTet(const G4String& pName,
           const G4ThreeVector& anchor,
           const G4ThreeVector& p1,
           const G4ThreeVector& p2,
           const G4ThreeVector& p3,
                 G4bool* degeneracyFlag = nullptr);

   ~G4UTet();

    void ComputeDimensions(      G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep);

    G4VSolid* Clone() const;

    inline G4GeometryType GetEntityType() const;

  public:   // without description

    G4UTet(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4UTet(const G4UTet& rhs);
    G4UTet& operator=(const G4UTet& rhs);
      // Copy constructor and assignment operator.

    void SetBoundingLimits(const G4ThreeVector& pMin, const G4ThreeVector& pMax);
    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const;

    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                                 G4double& pMin, G4double& pMax) const;

    G4Polyhedron* CreatePolyhedron() const;

    void SetVertices(const G4ThreeVector& anchor,
                     const G4ThreeVector& p1,
                     const G4ThreeVector& p2,
                     const G4ThreeVector& p3,
                     G4bool* degeneracyFlag = nullptr);
      // Set new position of the vertices.

    void GetVertices(G4ThreeVector& anchor,
                     G4ThreeVector& p1,
                     G4ThreeVector& p2,
                     G4ThreeVector& p3) const;
    std::vector<G4ThreeVector> GetVertices() const;
      // Return the four vertices of the shape.

    G4bool CheckDegeneracy(const G4ThreeVector& p0,
                           const G4ThreeVector& p1,
                           const G4ThreeVector& p2,
                           const G4ThreeVector& p3) const;
      // Return true if the tetrahedron is degenerate.

  private:

    G4ThreeVector fBmin, fBmax; // bounding box
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline G4GeometryType G4UTet::GetEntityType() const
{
  return "G4Tet";
}

#endif  // G4GEOM_USE_USOLIDS

#endif
