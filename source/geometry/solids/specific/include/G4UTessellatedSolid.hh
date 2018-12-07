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
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4UTesselladedSolid
//
// Class description:
//
//   Wrapper class for G4TessellatedSolid to make use of VecGeom TessellatedSolid.

// History:
// 11.01.18 G.Cosmo, CERN
// --------------------------------------------------------------------
#ifndef G4UTESSELLATEDSOLID_hh
#define G4UTESSELLATEDSOLID_hh

#include "G4UAdapter.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include <volumes/UnplacedTessellated.h>

#include "G4Polyhedron.hh"

class G4VFacet;

class G4UTessellatedSolid : public G4UAdapter<vecgeom::UnplacedTessellated>
{
  using Shape_t = vecgeom::UnplacedTessellated;
  using Base_t  = G4UAdapter<vecgeom::UnplacedTessellated>;

  public:  // with description

    G4UTessellatedSolid();
    G4UTessellatedSolid(const G4String& pName);
   ~G4UTessellatedSolid();

    G4bool AddFacet(G4VFacet* aFacet);
    G4VFacet* GetFacet(G4int i) const;

    G4int GetNumberOfFacets() const;

    inline G4GeometryType GetEntityType() const;

    void SetSolidClosed(const G4bool t);
    G4bool GetSolidClosed() const;

    void SetMaxVoxels(G4int);

    G4double GetMinXExtent() const;
    G4double GetMaxXExtent() const;
    G4double GetMinYExtent() const;
    G4double GetMaxYExtent() const;
    G4double GetMinZExtent() const;
    G4double GetMaxZExtent() const;
    G4int AllocatedMemoryWithoutVoxels();
    G4int AllocatedMemory();
    void DisplayAllocatedMemory();

  public:  // without description

    G4UTessellatedSolid(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4UTessellatedSolid( const G4UTessellatedSolid& source );
    G4UTessellatedSolid &operator=(const G4UTessellatedSolid& source);
      // Copy constructor and assignment operator.

    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const;
    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                                 G4double& pMin, G4double& pMax) const;  

    G4Polyhedron* CreatePolyhedron() const;

  private:

    std::vector<G4VFacet *>  fFacets;
    std::vector<G4ThreeVector>  fVertexList;
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline G4GeometryType G4UTessellatedSolid::GetEntityType() const
{
  return "G4TessellatedSolid";
}

#endif  // G4GEOM_USE_USOLIDS

#endif
