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
// G4UTesselladedSolid
//
// Class description:
//
// Wrapper class for G4TessellatedSolid to make use of VecGeom TessellatedSolid.

// Author: Gabriele Cosmo (CERN), 11.01.2018
// --------------------------------------------------------------------
#ifndef G4UTESSELLATEDSOLID_HH
#define G4UTESSELLATEDSOLID_HH

#include "G4UAdapter.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include <VecGeom/volumes/UnplacedTessellated.h>

#include "G4Polyhedron.hh"
#include "G4VFacet.hh"

/**
 * @brief G4UTessellatedSolid is a wrapper class for G4TessellatedSolid
 * to make use of VecGeom TessellatedSolid.
 */

class G4UTessellatedSolid : public G4UAdapter<vecgeom::UnplacedTessellated>
{
  using Shape_t = vecgeom::UnplacedTessellated;
  using Base_t  = G4UAdapter<vecgeom::UnplacedTessellated>;

  public:

    /**
     * Default Constructor.
     */
    G4UTessellatedSolid();

    /**
     * Constructor with solid's name.
     *  @param[in] name The name of the solid.
     */
    G4UTessellatedSolid(const G4String& pName);

    /**
     * Destructor. Clearing all allocated facets and data.
     */
   ~G4UTessellatedSolid() override;

    /**
     * Methods for adding or retrieving a facet given an index.
     */
    G4bool AddFacet(G4VFacet* aFacet);
    G4VFacet* GetFacet(G4int i) const;

    /**
     * Returns the total number of facets.
     */
    G4int GetNumberOfFacets() const;

    /**
     * Returns the type ID, "G4TessellatedSolid" of the solid.
     */
    inline G4GeometryType GetEntityType() const override;

    /**
     * Returns true as the solid has only planar faces.
     */
    inline G4bool IsFaceted() const override;

    /**
     * Modifier and accessor to close/finalise the solid.
     */
    void SetSolidClosed(const G4bool t);
    G4bool GetSolidClosed() const;

    /**
     * Allowing to tune the maximum number of voxels to use for optimisation.
     */
    void SetMaxVoxels(G4int);

    /**
     * Accessors.
     */
    G4double GetMinXExtent() const;
    G4double GetMaxXExtent() const;
    G4double GetMinYExtent() const;
    G4double GetMaxYExtent() const;
    G4double GetMinZExtent() const;
    G4double GetMaxZExtent() const;

    /**
     * Loggers reporting the total allocated memory.
     */
    G4int AllocatedMemoryWithoutVoxels();
    G4int AllocatedMemory();
    void DisplayAllocatedMemory();

    /**
     * Computes the bounding limits of the solid.
     *  @param[out] pMin The minimum bounding limit point.
     *  @param[out] pMax The maximum bounding limit point.
     */
    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const override;

    /**
     * Calculates the minimum and maximum extent of the solid, when under the
     * specified transform, and within the specified limits.
     *  @param[in] pAxis The axis along which compute the extent.
     *  @param[in] pVoxelLimit The limiting space dictated by voxels.
     *  @param[in] pTransform The internal transformation applied to the solid.
     *  @param[out] pMin The minimum extent value.
     *  @param[out] pMax The maximum extent value.
     *  @returns True if the solid is intersected by the extent region.
     */
    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                                 G4double& pMin, G4double& pMax) const override;  

    /**
     * Returns a generated polyhedron as graphical representations.
     */
    G4Polyhedron* CreatePolyhedron() const override;

    /**
     * Copy constructor and assignment operator.
     */
    G4UTessellatedSolid( const G4UTessellatedSolid& source );
    G4UTessellatedSolid& operator=(const G4UTessellatedSolid& source);

  private:

    std::vector<G4VFacet*> fFacets;
    std::vector<G4ThreeVector> fVertexList;
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline G4GeometryType G4UTessellatedSolid::GetEntityType() const
{
  return "G4TessellatedSolid";
}

inline G4bool G4UTessellatedSolid::IsFaceted() const
{
  return true;
}

#endif  // G4GEOM_USE_USOLIDS

#endif
