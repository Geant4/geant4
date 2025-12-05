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
// G4UPara
//
// Class description:
//
// Wrapper class for G4Para to make use of VecGeom Parallelepiped.

// Author: G.Cosmo (CERN), 13.09.2013
// --------------------------------------------------------------------
#ifndef G4UPARA_HH
#define G4UPARA_HH

#include "G4UAdapter.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include <VecGeom/volumes/UnplacedParallelepiped.h>

#include "G4Polyhedron.hh"

/**
 * @brief G4UPara is a wrapper class for G4Para to make use of
 * VecGeom Parallelepiped.
 */

class G4UPara : public G4UAdapter<vecgeom::UnplacedParallelepiped> 
{
  using Shape_t = vecgeom::UnplacedParallelepiped;
  using Base_t  = G4UAdapter<vecgeom::UnplacedParallelepiped>;

  public:

    /**
     * Constructs a parallelepiped, given a name and its parameters.
     *  @param[in] pName The name of the solid.
     *  @param[in] pDx Half-length in x.
     *  @param[in] pDy Half-length in y.
     *  @param[in] pDz Half-length in z.
     *  @param[in] pAlpha Angle formed by the Y axis and by the plane joining
     *             the centre of the faces parallel to the Z-X plane at -dy
     *             and +dy.
     *  @param[in] pTheta Polar angle of the line joining the centres of the
     *             faces at -dz and +dz in Z.
     *  @param[in] pPhi Azimuthal angle of the line joining the centres of
     *             the faces at -dz and +dz in Z.
     */
    G4UPara(const G4String& pName,
                  G4double pDx, G4double pDy, G4double pDz,
                  G4double pAlpha, G4double pTheta, G4double pPhi);

    /**
     * Constructs a parallelepiped, given a name and its 8 vertices.
     *  @param[in] pName The name of the solid.
     *  @param[in] pt Points of the 8 vertices.
     */
    G4UPara(const G4String& pName,
            const G4ThreeVector pt[8]);

    /**
     * Default destructor.
     */
    ~G4UPara() override = default;

    /**
     * Accessors.
     */
    G4double GetZHalfLength()  const;
    G4double GetYHalfLength()  const;
    G4double GetXHalfLength()  const;
    G4ThreeVector GetSymAxis() const;
    G4double GetTanAlpha()     const;

    /**
     * Accessors. Obtain (re)computed values of the original parameters.
     */
    G4double GetAlpha()  const;
    G4double GetTheta()  const;    
    G4double GetPhi()    const;
    // Obtain (re)computed values of original parameters
   
    /**
     * Modifiers.
     */
    void SetXHalfLength(G4double val);
    void SetYHalfLength(G4double val);
    void SetZHalfLength(G4double val);
    void SetAlpha(G4double alpha);
    void SetTanAlpha(G4double val);
    void SetThetaAndPhi(double pTheta, double pPhi);

    /**
     * Sets all parameters, as for constructor.
     */
    void SetAllParameters(G4double pDx, G4double pDy, G4double pDz,
                          G4double pAlpha, G4double pTheta, G4double pPhi);

    /**
     * Dispatch method for parameterisation replication mechanism and
     * dimension computation.
     */
    void ComputeDimensions(G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep) override;

    /**
     * Returns the type ID, "G4Para" of the solid.
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
     * Makes a clone of the object for use in multi-treading.
     *  @returns A pointer to the new cloned allocated solid.
     */
    G4VSolid* Clone() const override;

    /**
     * Returns a generated polyhedron as graphical representations.
     */
    G4Polyhedron* CreatePolyhedron() const override;

    /**
     * Copy constructor and assignment operator.
     */
    G4UPara(const G4UPara& rhs);
    G4UPara& operator=(const G4UPara& rhs);

  private:

    /**
     * Checks input parameters.
     */
    void CheckParameters();

    /**
     * Sets the side planes.
     */
    void MakePlanes();

  private:

    G4double fTalpha,fTthetaCphi,fTthetaSphi;
    struct { G4double a,b,c,d; } fPlanes[4];
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline G4GeometryType G4UPara::GetEntityType() const
{
  return "G4Para";
}

inline G4bool G4UPara::IsFaceted() const
{
  return true;
}

#endif  // G4GEOM_USE_USOLIDS

#endif
