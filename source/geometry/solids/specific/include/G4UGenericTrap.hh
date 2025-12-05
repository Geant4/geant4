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
// G4UGenericTrap
//
// Class description:
//
// Wrapper class for G4GenericTrap to make use of VecGeom GenericTrap.

// Author: Gabriele Cosmo (CERN), 30.10.2013
// --------------------------------------------------------------------
#ifndef G4UGENERICTRAP_HH
#define G4UGENERICTRAP_HH

#include "G4UAdapter.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include <VecGeom/volumes/UnplacedGenTrap.h>
#include "G4TwoVector.hh"

#include "G4Polyhedron.hh"

/**
 * @brief G4UGenericTrap is a wrapper class for G4GenericTrap
 * to make use of VecGeom GenericTrap.
 */

class G4UGenericTrap : public G4UAdapter<vecgeom::UnplacedGenTrap> 
{
  using Shape_t = vecgeom::UnplacedGenTrap;
  using Base_t  = G4UAdapter<vecgeom::UnplacedGenTrap>;

  public:

    /**
     * Constructs an generic trapezoid, given its vertices.
     *  @param[in] name The solid name.
     *  @param[in] halfZ Half length in Z.
     *  @param[in] vertices The (x,y) coordinates of the vertices.
     */
    G4UGenericTrap(const G4String& name, G4double halfZ,
                   const std::vector<G4TwoVector>& vertices);

    /**
     * Default destructor.
     */
    ~G4UGenericTrap() override = default;

    /**
     * Accessors.
     */
    G4double    GetZHalfLength() const;
    G4int       GetNofVertices() const;
    G4TwoVector GetVertex(G4int index) const;
    const std::vector<G4TwoVector>& GetVertices() const;
    G4double    GetTwistAngle(G4int index) const;
    G4bool      IsTwisted() const;
    G4int       GetVisSubdivisions() const;

    /**
     * Modifiers.
     */
    void        SetVisSubdivisions(G4int subdiv);
    void        SetZHalfLength(G4double);

    /**
     * Returns the type ID, "G4GenericTrap" of the solid.
     */
    inline G4GeometryType GetEntityType() const override;

    /**
     * Returns true as the solid has only planar faces.
     */
    inline G4bool IsFaceted() const override;

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
    G4UGenericTrap( const G4UGenericTrap& source );
    G4UGenericTrap& operator=(const G4UGenericTrap& source);

  private:

    /**
     * Initialises data. Used in constructor.
     */
    void Initialise(const std::vector<G4TwoVector>& v);

  private:

    G4int fVisSubdivisions;
    std::vector<G4TwoVector> fVertices;

};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline G4GeometryType G4UGenericTrap::GetEntityType() const
{
  return "G4GenericTrap";
}

inline G4bool G4UGenericTrap::IsFaceted() const
{
  return true;
}

#endif  // G4GEOM_USE_USOLIDS

#endif
