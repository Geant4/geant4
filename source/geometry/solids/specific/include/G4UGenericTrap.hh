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
// G4UGenericTrap
//
// Class description:
//
//   Wrapper class for G4UGenericTrap to make use of it from USolids module.

// History:
// 30.10.13 G.Cosmo, CERN/PH
// --------------------------------------------------------------------
#ifndef G4UGENERICTRAP_hh
#define G4UGENERICTRAP_hh

#include "G4UAdapter.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include <volumes/UnplacedGenTrap.h>
#include "G4TwoVector.hh"

#include "G4Polyhedron.hh"

class G4UGenericTrap : public G4UAdapter<vecgeom::UnplacedGenTrap> 
{
  using Shape_t = vecgeom::UnplacedGenTrap;
  using Base_t  = G4UAdapter<vecgeom::UnplacedGenTrap>;

  public:  // with description

    G4UGenericTrap(const G4String& name, G4double halfZ,
                   const std::vector<G4TwoVector>& vertices);

   ~G4UGenericTrap();

    G4double    GetZHalfLength() const;
    G4int       GetNofVertices() const;
    G4TwoVector GetVertex(G4int index) const;
    const std::vector<G4TwoVector>& GetVertices() const;
    G4double    GetTwistAngle(G4int index) const;
    G4bool      IsTwisted() const;
    G4int       GetVisSubdivisions() const;
    void        SetVisSubdivisions(G4int subdiv);
    void        SetZHalfLength(G4double);
    void Initialise(const std::vector<G4TwoVector>& v);

    inline G4GeometryType GetEntityType() const;

  public:  // without description

    G4UGenericTrap(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4UGenericTrap( const G4UGenericTrap& source );
    G4UGenericTrap &operator=(const G4UGenericTrap& source);
      // Copy constructor and assignment operator.

    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const;

    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                                 G4double& pMin, G4double& pMax) const;

    G4Polyhedron* CreatePolyhedron() const;

  private:

    G4int                    fVisSubdivisions;
    std::vector<G4TwoVector> fVertices;

};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline G4GeometryType G4UGenericTrap::GetEntityType() const
{
  return "G4GenericTrap";
}

#endif  // G4GEOM_USE_USOLIDS

#endif
