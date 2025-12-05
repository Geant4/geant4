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
// G4SubtractionSolid
//
// Class description:
//
// Class for description of subtraction of two solids: A - B.

// Author: Vladimir Grichine (CERN), 14.10.1998 - First implementation
// --------------------------------------------------------------------
#ifndef G4SUBTRACTIONSOLID_HH
#define G4SUBTRACTIONSOLID_HH

#include "G4BooleanSolid.hh"
#include "G4VSolid.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4AffineTransform.hh"

/**
 * @brief G4SubtractionSolid is a solid describing the Boolean subtraction
 * of two solids.
 */

class G4SubtractionSolid : public G4BooleanSolid
{
  public:

    /**
     * Constructor of a Boolean subtraction between two solids with no
     * displacement.
     *  @param[in] pName The name of the Boolean composition.
     *  @param[in] pSolidA Pointer to the first reference solid.
     *  @param[in] pSolidB Pointer to the second solid to form the composition.
     */
    G4SubtractionSolid(  const G4String& pName,
                               G4VSolid* pSolidA ,
                               G4VSolid* pSolidB   ) ;

    /**
     * Constructor of a Boolean subtraction between two solids with rotation
     * and translation, used to transform the coordinate system of the second
     * solid to the coordinate system of the first solid.
     *  @param[in] pName The name of the Boolean composition.
     *  @param[in] pSolidA Pointer to the first reference solid.
     *  @param[in] pSolidB Pointer to the second solid to form the composition.
     *  @param[in] rotMatrix Pointer to the rotation vector.
     *  @param[in] transVector The translation vector.
     */
    G4SubtractionSolid(  const G4String& pName,
                               G4VSolid* pSolidA ,
                               G4VSolid* pSolidB ,
                               G4RotationMatrix* rotMatrix,
                         const G4ThreeVector& transVector   ) ;

    /**
     * Constructor of a Boolean subtraction between two solids with a
     * transformation that moves the second solid from its desired position
     * to its standard position.
     *  @param[in] pName The name of the Boolean composition.
     *  @param[in] pSolidA Pointer to the first reference solid.
     *  @param[in] pSolidB Pointer to the second solid to form the composition.
     *  @param[in] transform The composed 3D transformation.
     */
    G4SubtractionSolid(  const G4String& pName,
                               G4VSolid* pSolidA ,
                               G4VSolid* pSolidB ,
                         const G4Transform3D& transform   ) ;

    /**
     * Default destructor.
     */
    ~G4SubtractionSolid() override = default ;

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4SubtractionSolid(__void__&);

    /**
     * Copy constructor and assignment operator.
     */
    G4SubtractionSolid(const G4SubtractionSolid& rhs);
    G4SubtractionSolid& operator=(const G4SubtractionSolid& rhs);

    /**
     * Returns the type ID, "G4SubtractionSolid" of the solid.
     */
    G4GeometryType GetEntityType() const override ;

    /**
     * Makes a clone of the object for use in multi-treading.
     *  @returns A pointer to the new cloned allocated solid.
     */
    G4VSolid* Clone() const override;

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
                                  G4double& pMin, G4double& pMax) const override ;
       
    /**
     * Returns if the given point "p" is inside or not the solid.
     */
    EInside Inside( const G4ThreeVector& p ) const override ;

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
     * Methods for creating graphical representations (i.e. for visualisation).
     */
    void DescribeYourselfTo ( G4VGraphicsScene& scene ) const override ;
    G4Polyhedron* CreatePolyhedron () const override ;

    /**
     * Returns an estimate of the capacity of the Boolean composition.
     */
    G4double GetCubicVolume() final;
};

#endif
