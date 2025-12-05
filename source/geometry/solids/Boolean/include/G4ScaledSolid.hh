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
// G4ScaledSolid
//
// Class description:
//
// A scaled solid is a solid that has been scaled in dimensions
// in X, Y or Z, from its original description.

// Author: Gabriele Cosmo (CERN), 27.10.2015 - Created
// --------------------------------------------------------------------
#ifndef G4SCALEDSOLID_HH
#define G4SCALEDSOLID_HH

#include "G4VSolid.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4AffineTransform.hh"

class G4ScaleTransform;

/**
 * @brief G4ScaledSolid is a solid that has been scaled in dimensions
 * in X, Y or Z, from its original description.
 */

class G4ScaledSolid : public G4VSolid
{
  public:

    /**
     * Constructor of a solid with scaled transformation.
     *  @param[in] pName The name of the solid.
     *  @param[in] pSolid Pointer to the original reference solid.
     *  @param[in] pScale The scaling transformation.
     */
    G4ScaledSolid( const G4String& pName,
                         G4VSolid* pSolid ,
                   const G4Scale3D& pScale  );

    /**
     * The destructor, clearing the cached transformation.
     */
    ~G4ScaledSolid() override;

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4ScaledSolid(__void__&);

    /**
     * Copy constructor and assignment operator.
     */
    G4ScaledSolid(const G4ScaledSolid& rhs);
    G4ScaledSolid& operator=(const G4ScaledSolid& rhs);

    /**
     * Returns if the given point "p" is inside or not the solid.
     */
    EInside Inside( const G4ThreeVector& p ) const override;

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
     * Returns the outwards pointing unit normal of the shape for the
     * surface closest to the point at offset "p".
     */
    G4ThreeVector SurfaceNormal( const G4ThreeVector& p ) const override;

    /**
     * Returns the distance along the normalised vector "v" to the shape,
     * from the point at offset "p". If there is no intersection, return
     * kInfinity. The first intersection resulting from leaving a
     * surface/volume is discarded. Hence, it is tolerant of points on
     * the surface of the shape.
     */
    G4double DistanceToIn( const G4ThreeVector& p,
                           const G4ThreeVector& v  ) const override;

    /**
     * Calculates the safety distance to the nearest surface of a shape from
     * an outside point. The distance can be an underestimate.
     */
    G4double DistanceToIn( const G4ThreeVector& p) const override;

    /**
     * Returns the distance along the normalised vector "v" to the shape,
     * from a point at an offset "p" inside or on the surface of the shape.
     * Intersections with surfaces, when the point is < Tolerance/2 from a
     * surface must be ignored. Must be called as solid.DistanceToOut(p,v)
     * or by specifying all the parameters.
     *  @param[in] p The reference point in space.
     *  @param[in] v The normalised direction.
     *  @param[in] calcNorm Flag to enable the normal computation or not.
     *  @param[out] validNorm Set to true if the solid lies entirely behind
     *              or on the exiting surface (calcNorm must be true, otherwise
     *              it is unused).
     *  @param[out] n The exiting outwards normal vector (undefined Magnitude).
     *              (calcNorm must be true, otherwise it is unused). 
     *  @returns The distance value to exit the volume.
     */
    G4double DistanceToOut( const G4ThreeVector& p,
                            const G4ThreeVector& v,
                            const G4bool calcNorm = false,
                                  G4bool* validNorm = nullptr,
                                  G4ThreeVector* n = nullptr ) const override;

    /**
     * Calculates the safety distance to the nearest surface of a shape from
     * an inside point "p". The distance can be an underestimate.
     */
    G4double DistanceToOut( const G4ThreeVector& p ) const override;

    /**
     * Throws an exception as paramterisations are not allowed for these solids.
     */
    void ComputeDimensions(       G4VPVParameterisation* p,
                            const G4int n,
                            const G4VPhysicalVolume* pRep ) override;

    /**
     * Methods returning an estimation of the solid volume (capacity) and
     * surface area, in internal units.
     */
    G4double GetCubicVolume() override;
    G4double GetSurfaceArea() override;

    /**
     * Returns a random point located on the surface of the solid.
     * Points returned may not necessarily be uniformly distributed.
     */
    G4ThreeVector GetPointOnSurface() const override;

    /**
     * Returns the number of constituents of the solid.
     * For non-Boolean solids the return value is one.
     */
    G4int GetNumOfConstituents() const override;

    /**
     * Returns true if the solid has only planar faces, false otherwise.
     */
    G4bool IsFaceted() const override;

    /**
     * Accessor and setter for the scaling transformation.
     */
    G4Scale3D GetScaleTransform() const; 
    void SetScaleTransform(const G4Scale3D& scale); 

    /**
     * Returns a pointer to the original not scaled  solid.
     */
    G4VSolid* GetUnscaledSolid() const;

    /**
     * Returns the type ID, "G4ScaledSolid" of the solid.
     */
    G4GeometryType  GetEntityType() const override;

    /**
     * Makes a clone of the object for use in multi-treading.
     *  @returns A pointer to the new cloned allocated solid.
     */
    G4VSolid* Clone() const override;

    /**
     * Streams the object contents to an output stream.
     */
    std::ostream& StreamInfo(std::ostream& os) const override;

    /**
     * Methods for creating graphical representations (i.e. for visualisation).
     */
    void DescribeYourselfTo ( G4VGraphicsScene& scene ) const override;
    G4Polyhedron* CreatePolyhedron () const override;
    G4Polyhedron* GetPolyhedron    () const override;

  private:

    G4VSolid* fPtrSolid = nullptr;
    G4ScaleTransform* fScale = nullptr;
    G4double fCubicVolume = -1.0;
    G4double fSurfaceArea = -1.0;
    mutable G4bool fRebuildPolyhedron = false;
    mutable G4Polyhedron* fpPolyhedron = nullptr;
};

#endif
