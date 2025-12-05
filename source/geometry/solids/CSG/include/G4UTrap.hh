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
// G4UTrap
//
// Class description:
//
// Wrapper class for G4Trap to make use of VecGeom Trapezoid.

// Author: G.Cosmo (CERN), 13.09.2013
// --------------------------------------------------------------------
#ifndef G4UTRAP_HH
#define G4UTRAP_HH

#include "G4UAdapter.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include <VecGeom/volumes/UnplacedTrapezoid.h>

#include "G4Polyhedron.hh"

/**
 * @brief G4UTrap is a wrapper class for G4Trap to make use of VecGeom Trapezoid.
 */

class G4UTrap : public G4UAdapter<vecgeom::UnplacedTrapezoid>
{
  using Shape_t = vecgeom::UnplacedTrapezoid;
  using Base_t = G4UAdapter<vecgeom::UnplacedTrapezoid>;

  public:

    /**
     * The most general constructor for G4Trap which prepares plane
     * equations and corner coordinates from parameters.
     *  @param[in] pName The name of the solid.
     *  @param[in] pDz Half-length along the Z-axis.
     *  @param[in] pTheta Polar angle of the line joining the centres
     *             of the faces at -/+pDz.
     *  @param[in] pPhi Azimuthal angle of the line joining the centre
     *             of the face at -pDz to the centre of the face at +pDz.
     *  @param[in] pDy1 Half-length along Y of the face at -pDz.
     *  @param[in] pDx1 Half-length along X of the side at y=-pDy1
     *             of the face at -pDz.
     *  @param[in] pDx2 Half-length along X of the side at y=+pDy1
     *             of the face at -pDz.
     *  @param[in] pAlp1 Angle with respect to the Y axis from the centre of the
     *             side at y=-pDy1 to the centre at y=+pDy1 of the face at -pDz.
     *  @param[in] pDy2 Half-length along Y of the face at +pDz.
     *  @param[in] pDx3 Half-length along X of the side at y=-pDy2
     *             of the face at +pDz.
     *  @param[in] pDx4 Half-length along X of the side at y=+pDy2
     *             of the face at +pDz.
     *  @param[in] pAlp2 Angle with respect to the Y axis from the centre of the
     *             side at y=-pDy2 to the centre at y=+pDy2 of the face at +pDz.
     */
    G4UTrap( const G4String& pName,
                   G4double pDz,
                   G4double pTheta, G4double pPhi,
                   G4double pDy1, G4double pDx1, G4double pDx2,
                   G4double pAlp1,
                   G4double pDy2, G4double pDx3, G4double pDx4,
                   G4double pAlp2 );

    /**
     * Prepares plane equations and parameters from corner coordinates.
     *  @param[in] pName The name of the solid.
     *  @param[in] pt Points of the 8 vertices.
     */
    G4UTrap( const G4String& pName,
             const G4ThreeVector pt[8] ) ;

    /**
     * Constructor for Right Angular Wedge from STEP (assumes pLTX<=pX).
     *  @param[in] pName The name of the solid.
     *  @param[in] pZ Length along Z.
     *  @param[in] pY Length along Y.
     *  @param[in] pX Length along X at the wider side.
     *  @param[in] pLTX Length along X at the narrower side (plTX<=pX).
     */
    G4UTrap( const G4String& pName,
                   G4double pZ,
                   G4double pY,
                   G4double pX, G4double pLTX );

    /**
     * Constructor for G4Trd.
     *  @param[in] pName The name of the solid.
     *  @param[in] pDx1 Half-length along X at the surface positioned at -dz.
     *  @param[in] pDx2 Half-length along X at the surface positioned at +dz.
     *  @param[in] pDy1 Half-length along Y at the surface positioned at -dz.
     *  @param[in] pDy2 Half-length along Y at the surface positioned at +dz.
     *  @param[in] pDz Half-length along Z axis.
     */
    G4UTrap( const G4String& pName,
                   G4double pDx1,  G4double pDx2,
                   G4double pDy1,  G4double pDy2,
                   G4double pDz );

    /**
     * Constructor for G4Para.
     *  @param[in] pName The name of the solid.
     *  @param[in] pDx Half-length in X.
     *  @param[in] pDy Half-length in Y.
     *  @param[in] pDz Half-length in Z.
     *  @param[in] pAlpha Angle formed by the Y axis and the plane joining the
     *             centre of the faces parallel to the Z-X plane at -dy and +dy.
     *  @param[in] pTheta Polar angle of the line joining the centres of the
     *             faces at -dz and +dz in Z.
     *  @param[in] pPhi Azimuthal angle of the line joining the centres of
     *             the faces at -dz and +dz in Z.
     */
    G4UTrap(const G4String& pName,
                  G4double pDx, G4double pDy, G4double pDz,
                  G4double pAlpha, G4double pTheta, G4double pPhi );

    /**
     * Constructor for "nominal" G4Trap whose parameters are to be set
     * by a G4VPVParamaterisation later on.
     *  @param[in] pName The name of the solid.
     */
    G4UTrap( const G4String& pName );

    /**
     * Default destructor.
     */
    ~G4UTrap() override = default;

    /**
     * Dispatch method for parameterisation replication mechanism and
     * dimension computation.
     */
    void ComputeDimensions(G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep) override;

    /**
     * Makes a clone of the object for use in multi-treading.
     *  @returns A pointer to the new cloned allocated solid.
     */
    G4VSolid* Clone() const override;

    using Base_t::GetTanAlpha1;
    using Base_t::GetTanAlpha2;

    /**
     * Accessors. Returning the coordinates of a unit vector along a straight
     * line joining centers of -/+fDz planes.
     */
    G4double GetZHalfLength()  const;
    G4double GetYHalfLength1() const;
    G4double GetXHalfLength1() const;
    G4double GetXHalfLength2() const;
    G4double GetTanAlpha1()    const;
    G4double GetYHalfLength2() const;
    G4double GetXHalfLength3() const;
    G4double GetXHalfLength4() const;
    G4double GetTanAlpha2()    const;

    /**
     * More accessors.
     */
    TrapSidePlane GetSidePlane(G4int n) const;
    G4ThreeVector GetSymAxis() const;

    /**
     * Accessors obtaining (re)computed values of the original parameters.
     */
    G4double GetPhi()    const;
    G4double GetTheta()  const;
    G4double GetAlpha1() const;
    G4double GetAlpha2() const;

    /**
     * Sets all parameters, as for constructor. Checks and sets half-widths
     * as well as angles. Makes a final check of co-planarity.
     */
    void SetAllParameters(G4double pDz, G4double pTheta, G4double pPhi,
                          G4double pDy1, G4double pDx1, G4double pDx2,
                          G4double pAlp1,
                          G4double pDy2, G4double pDx3, G4double pDx4,
                          G4double pAlp2);

    /**
     * Returns the type ID, "G4Trap" of the solid.
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
    G4UTrap(const G4UTrap& rhs);
    G4UTrap& operator=(const G4UTrap& rhs);

  private:

    /**
     * Sets parameters using eight vertices.
     */
    void SetPlanes(const G4ThreeVector pt[8]);

    /**
     * Checks dimensions.
     */
    void CheckParameters() const;

    /**
     * Computes coordinates of vertices.
     */
    void GetVertices(G4ThreeVector pt[8]) const;

    /**
     * Checks planarity of lateral planes.
     */
    void CheckPlanarity(const G4ThreeVector pt[8]) const;
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline G4GeometryType G4UTrap::GetEntityType() const
{
  return "G4Trap";
}

inline G4bool G4UTrap::IsFaceted() const
{
  return true;
}

#endif  // G4GEOM_USE_USOLIDS

#endif
