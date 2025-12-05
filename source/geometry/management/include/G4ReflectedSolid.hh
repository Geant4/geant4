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
// G4ReflectedSolid
//
// Class description:
//
// A Reflected solid is a solid that has been shifted from its original
// frame of reference to a new reflected one.

// Author: Vladimir Grichine (CERN), 23.07.2001 - Created
// --------------------------------------------------------------------
#ifndef G4ReflectedSolid_hh
#define G4ReflectedSolid_hh

#include "G4VSolid.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"

/**
 * @brief G4ReflectedSolid is a solid that has been shifted from its original
 * frame of reference to a new reflected one.
 */

class G4ReflectedSolid : public G4VSolid
{
  public:

    /**
     * Constructor for G4ReflectedSolid. For use in instantiating
     * a transient instance.
     *  @param[in] pName The solid's name.
     *  @param[in] pSolid The original primitive being reflected.
     *  @param[in] transform The associated transformation.
     */
    G4ReflectedSolid( const G4String& pName,
                            G4VSolid* pSolid ,
                      const G4Transform3D& transform ) ;

    /**
     * Destructor.
     */
    ~G4ReflectedSolid() override;

    // Includes all the methods that a solid requires.

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
    G4bool CalculateExtent( const EAxis pAxis,
                            const G4VoxelLimits& pVoxelLimit,
                            const G4AffineTransform& pTransform,
                                  G4double& pMin, G4double& pMax) const override;

    /**
     * Concrete implementations of the expected query interfaces for
     * solids, as defined in the base class G4VSolid.
     */
    EInside Inside( const G4ThreeVector& p ) const override; 
    G4ThreeVector SurfaceNormal( const G4ThreeVector& p ) const override;
    G4double DistanceToIn( const G4ThreeVector& p,
                           const G4ThreeVector& v ) const override;
    G4double DistanceToIn( const G4ThreeVector& p) const override;
    G4double DistanceToOut( const G4ThreeVector& p,
                            const G4ThreeVector& v,
                            const G4bool calcNorm = false,
                                  G4bool* validNorm = nullptr,
                                  G4ThreeVector* n = nullptr ) const override;
    G4double DistanceToOut( const G4ThreeVector& p ) const override;

    /**
     * Dispatch method for parameterisation replication mechanism and
     * dimension computation.
     */
    void ComputeDimensions( G4VPVParameterisation* p,
                            const G4int n,
                            const G4VPhysicalVolume* pRep ) override;

    /**
     * Returning an estimation of the solid volume (capacity) and
     * surface area, in internal units.
     */
    G4double GetCubicVolume() override;
    G4double GetSurfaceArea() override;

    /**
     * Returns a random point located and uniformly distributed on the
     * surface of the solid.
     */
    G4ThreeVector GetPointOnSurface() const override;

    /**
     * Returns the number of constituent solids (in case Boolean).
     */
    G4int GetNumOfConstituents() const override;

    /**
     * Returns true as the solid has only planar faces.
     */
    G4bool IsFaceted() const override;

    /**
     * Makes a clone of the object for use in multi-treading.
     *  @returns A pointer to the new cloned allocated solid.
     */
    G4VSolid* Clone() const override;

    /**
     * Returns a random point located and uniformly distributed on the
     * surface of the solid.
     */
    G4GeometryType  GetEntityType() const override;

    /**
     * If the Solid is a G4ReflectedSolid, return a self pointer else
     * return nullptr.
     */
    virtual const G4ReflectedSolid* GetReflectedSolidPtr() const;
    virtual G4ReflectedSolid* GetReflectedSolidPtr();

    /**
     * Returns a pointer to the original solid primitive.
     */
    G4VSolid* GetConstituentMovedSolid() const;

    /**
     * Accessors and modifier for the transformation.
     */
    G4Transform3D GetTransform3D() const; 
    G4Transform3D GetDirectTransform3D() const; 
    void SetDirectTransform3D(G4Transform3D&);

    /**
     * Streams the object contents to an output stream.
     */
    std::ostream& StreamInfo(std::ostream& os) const override;

    /**
     * Copy constructor and assignment operator.
     */
    G4ReflectedSolid(const G4ReflectedSolid& rhs);
    G4ReflectedSolid& operator=(const G4ReflectedSolid& rhs);

    /**
     * Methods for creating graphical representations (i.e. for visualisation).
     */
    void DescribeYourselfTo(G4VGraphicsScene& scene) const override;
    G4Polyhedron* CreatePolyhedron() const override;
    G4Polyhedron* GetPolyhedron() const override;

  protected:

    G4VSolid*          fPtrSolid = nullptr;
    G4Transform3D*     fDirectTransform3D = nullptr;

    /** Caches for the reflected G4Polyhedron. */
    mutable G4bool fRebuildPolyhedron = false;
    mutable G4Polyhedron* fpPolyhedron = nullptr;
};

#endif
