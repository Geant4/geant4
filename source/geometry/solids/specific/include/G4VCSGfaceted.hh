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
// G4VCSGfaceted
//
// Class description:
//
// Virtual class defining CSG-like type shape that is built entirely
// of G4CSGface faces.

// Author: David C. Williams (UCSC), 1998 - Created
// --------------------------------------------------------------------
#ifndef G4VCSGFACETED_HH
#define G4VCSGFACETED_HH 1

#include "G4VSolid.hh"

class G4VCSGface;
class G4VisExtent;

/**
 * @brief G4VCSGfaceted is a virtual class defining a CSG-like type shape
 * that is built entirely of G4CSGface faces.
 */

class G4VCSGfaceted : public G4VSolid 
{
  public:

    /**
     * Constructor taking a 'name'.
     */
    G4VCSGfaceted( const G4String& name );

    /**
     * Destructor.
     */
    ~G4VCSGfaceted() override;
  
    /**
     * Copy constructor and assignment operator.
     */
    G4VCSGfaceted( const G4VCSGfaceted& source );
    G4VCSGfaceted& operator=( const G4VCSGfaceted& source );
  
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
    G4bool CalculateExtent( const EAxis pAxis,
                            const G4VoxelLimits& pVoxelLimit,
                            const G4AffineTransform& pTransform,
                                  G4double& pmin, G4double& pmax) const override;
  
    /**
     * Concrete implementations of the expected query interfaces for
     * solids, as defined in G4VSolid.
     */
    EInside Inside( const G4ThreeVector& p ) const override;
    G4ThreeVector SurfaceNormal( const G4ThreeVector& p ) const override;
    G4double DistanceToIn( const G4ThreeVector& p,
                           const G4ThreeVector& v ) const override;
    G4double DistanceToIn( const G4ThreeVector& p ) const override;
    G4double DistanceToOut( const G4ThreeVector& p,
                            const G4ThreeVector& v,
                            const G4bool calcNorm = false,
                                  G4bool* validNorm = nullptr,
                                  G4ThreeVector* n = nullptr ) const override;
    G4double DistanceToOut( const G4ThreeVector& p ) const override;

    /**
     * Returns the type ID, "G4CSGfaceted" of the solid.
     */
    G4GeometryType GetEntityType() const override;

    /**
     * Streams the object contents to an output stream.
     */
    std::ostream& StreamInfo(std::ostream& os) const override;

    /**
     * Returns a pointer to a generated polyhedron used for visualisation.
     */
    G4Polyhedron* CreatePolyhedron() const override = 0;

    /**
     * Methods for creating graphical representations (i.e. for visualisation).
     */
    void DescribeYourselfTo( G4VGraphicsScene& scene ) const override;
    G4VisExtent GetExtent() const override;
    G4Polyhedron* GetPolyhedron () const override;

    /**
     * Accessors and modifiers for capacity and area computation.
     */
    G4int GetCubVolStatistics() const;
    G4double GetCubVolEpsilon() const;
    void SetCubVolStatistics(G4int st);
    void SetCubVolEpsilon(G4double ep);
    G4int GetAreaStatistics() const;
    G4double GetAreaAccuracy() const;
    void SetAreaStatistics(G4int st);
    void SetAreaAccuracy(G4double ep);

    /**
     * Returning an estimation of the solid volume (capacity) and
     * surface area, in internal units. Caches the computed value
     * once computed the first time.
     */
    G4double GetCubicVolume() override;
    G4double GetSurfaceArea() override;

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4VCSGfaceted(__void__&);

  protected:

    /**
     * Protected method used in DistanceToIn() and DistanceToOut().
     */
    virtual G4double DistanceTo( const G4ThreeVector& p,
                                 const G4bool outgoing ) const;

    /**
     * Returns a random point located on the surface of the solid 
     * in case of generic Polycone or generic Polyhedra.
     */
    G4ThreeVector GetPointOnSurfaceGeneric()const;

    /**
     * Copy parameters from other solid or reset them.
     * Used in copy constructor and assignment operator.
     */
    void CopyStuff( const G4VCSGfaceted& source );
    void DeleteStuff();

  protected:

    G4int numFace = 0;
    G4VCSGface **faces = nullptr;
    G4double fCubicVolume = 0.0;
    G4double fSurfaceArea = 0.0;
    mutable G4bool fRebuildPolyhedron = false;
    mutable G4Polyhedron* fpPolyhedron = nullptr;

  private:

    /** Statistics, error accuracy for volume estimation. */
    G4int fStatistics;
    G4double fCubVolEpsilon;
    G4double fAreaAccuracy;
};

#endif
