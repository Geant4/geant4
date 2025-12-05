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
// G4DisplacedSolid
//
// Class description:
//
// A displaced solid is a solid that has been shifted from its original
// frame of reference to a new one. It is meant to be used **internally only**
// for simplifying the implementation of "Boolean solids". 

// Author: Vladimir Grichine (CERN), 28.10.1998 - Created.
// --------------------------------------------------------------------
#ifndef G4DISPLACEDSOLID_HH
#define G4DISPLACEDSOLID_HH

#include "G4VSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4AffineTransform.hh"

/**
 * @brief G4DisplacedSolid is a solid that has been shifted from its original
 * frame of reference to a new one. It is meant to be used **internally only**,
 * for simplifying the implementation of "Boolean solids".
 */

class G4DisplacedSolid : public G4VSolid
{
  public:

    /**
     * Constructor of a displaced solid rotation and translation vectors.
     *  @param[in] pName The name of the diplaced solid.
     *  @param[in] pSolid Pointer to the original reference solid.
     *  @param[in] rotMatrix Pointer to the rotation vector.
     *  @param[in] transVector The translation vector.
     */
    G4DisplacedSolid( const G4String& pName,
                            G4VSolid* pSolid ,
                            G4RotationMatrix* rotMatrix,
                      const G4ThreeVector& transVector  ) ;

    /**
     * Constructor of a displaced solid with a transformation.
     *  @param[in] pName The name of the displaced solid.
     *  @param[in] pSolid Pointer to the original reference solid.
     *  @param[in] transform The composed 3D transformation.
     */
    G4DisplacedSolid( const G4String& pName,
                            G4VSolid* pSolid ,
                      const G4Transform3D& transform  ) ;

    /**
     * Constructor for use in instantiating a transient instance from a
     * persistent one.
     *  @param[in] pName The name of the displaced solid.
     *  @param[in] pSolid Pointer to the original reference solid.
     *  @param[in] directTransform The internal transformation.
     */
    G4DisplacedSolid( const G4String& pName,
                            G4VSolid* pSolid ,
                      const G4AffineTransform directTransform );

    /**
     * Destructor. Deletes all cached transformations.
     */
    ~G4DisplacedSolid() override ;

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4DisplacedSolid(__void__&);

    /**
     * Copy constructor and assignment operator.
     */
    G4DisplacedSolid(const G4DisplacedSolid& rhs);
    G4DisplacedSolid& operator=(const G4DisplacedSolid& rhs);

    /**
     * Returns if the given point "p" is inside or not the solid.
     */
    EInside Inside( const G4ThreeVector& p ) const override ; 

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
                                 G4double& pMin, G4double& pMax) const override ;

    /**
     * Returns the outwards pointing unit normal of the shape for the
     * surface closest to the point at offset "p".
     */
    G4ThreeVector SurfaceNormal( const G4ThreeVector& p ) const override ;

    /**
     * Returns the distance along the normalised vector "v" to the shape,
     * from the point at offset "p". If there is no intersection, return
     * kInfinity. The first intersection resulting from leaving a
     * surface/volume is discarded. Hence, it is tolerant of points on
     * the surface of the shape.
     */
    G4double DistanceToIn( const G4ThreeVector& p,
                           const G4ThreeVector& v  ) const override ;

    /**
     * Calculates the safety distance to the nearest surface of a shape from
     * an outside point. The distance can be an underestimate.
     */
    G4double DistanceToIn( const G4ThreeVector& p) const override ;

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
                                  G4ThreeVector* n = nullptr ) const override ;

    /**
     * Calculates the safety distance to the nearest surface of a shape from
     * an inside point "p". The distance can be an underestimate.
     */
    G4double DistanceToOut( const G4ThreeVector& p ) const override ;


    /**
     * Throws an exception as paramterisations are not allowed for these solids.
     */
    void ComputeDimensions(       G4VPVParameterisation* p,
                            const G4int n,
                            const G4VPhysicalVolume* pRep ) override ;

    /**
     * Deletes cached transformations. Used in destructor.
     */
    void CleanTransformations();

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
     * Returns the type ID, "G4DisplacedSolid" of the solid.
     */
    G4GeometryType  GetEntityType() const override;

    /**
     * Makes a clone of the object for use in multi-treading.
     *  @returns A pointer to the new cloned allocated solid.
     */
    G4VSolid* Clone() const override;

    /**
     * If the Solid is a "G4DisplacedSolid", return a self pointer else
     * return nullptr.
     */
    const G4DisplacedSolid* GetDisplacedSolidPtr() const override;
          G4DisplacedSolid* GetDisplacedSolidPtr() override;

    /**
     * Returns a pointer to the original not displaced solid.
     */
    G4VSolid* GetConstituentMovedSolid() const;

    /**
     * Accessor/modifier for the associated internal transformation.
     */
    G4AffineTransform GetTransform() const; 
    void SetTransform(G4AffineTransform& ); 

    /**
     * Accessor/modifier for the associated internal transformation, as above.
     */
    G4AffineTransform GetDirectTransform() const; 
    void SetDirectTransform(G4AffineTransform&); 

    /**
     * Get/Set the rotation/translation, as applied to the frame of reference.
     */
    G4RotationMatrix GetFrameRotation() const;
    void SetFrameRotation(const G4RotationMatrix&);
    G4ThreeVector GetFrameTranslation() const; 
    void SetFrameTranslation(const G4ThreeVector&); 

    /**
     * Get/Set the rotation/translation, as applied to the object.
     */
    G4RotationMatrix GetObjectRotation() const;
    void SetObjectRotation(const G4RotationMatrix&);
    G4ThreeVector GetObjectTranslation() const; 
    void SetObjectTranslation(const G4ThreeVector&); 

    /**
     * Streams the object contents to an output stream.
     */
    std::ostream& StreamInfo(std::ostream& os) const override;

    /**
     * Methods for creating graphical representations (i.e. for visualisation).
     */
    void DescribeYourselfTo ( G4VGraphicsScene& scene ) const override ;
    G4Polyhedron* CreatePolyhedron () const override ;
    G4Polyhedron* GetPolyhedron () const override ;

  protected:

    G4VSolid* fPtrSolid = nullptr;
    G4AffineTransform* fPtrTransform = nullptr;
    G4AffineTransform* fDirectTransform = nullptr;
    mutable G4bool fRebuildPolyhedron = false;
    mutable G4Polyhedron* fpPolyhedron = nullptr;
};

#endif
