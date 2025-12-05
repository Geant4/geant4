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
// G4VCSGface
//
// Class description:
//
// Definition of the virtual base class G4VCSGface, one side (or face)
// of a CSG-like solid. It should be possible to build a CSG entirely
// out of connecting CSG faces.
//
// Each face has an inside and outside surface, the former represents
// the inside of the volume, the latter, the outside.
//
// -------------------------------------------------------------------
//
// Implementation notes:
//   * distance.
//        The meaning of distance includes the boundaries of the face.
//        For example, for a rectangular, planer face:
//
//               A   |  B           | C
//                   |              |
//              -------+--------------+-----
//               D   |  I           | E
//                   |              |
//              -------+--------------+-----
//               F   |  G           | H
//                   |              |
//       
//        A, C, F, and H: closest distance is the distance to
//        the adjacent corner.
//
//        B, D, E, and G: closest distance is the distance to
//        the adjacent line.
//
//        I: normal distance to plane
//
//        For non-planer faces, one can use the normal to decide when
//        a point falls off the edge and then act accordingly.
//
//
// Usage:
//
//   A CSG shape can be defined by putting together any number of generic
//   faces, as long as the faces cover the entire surface of the shape
//   without overlapping.
//
//   G4VSolid::CalculateExtent
//
//   Define unit vectors along the specified transform axis.
//   Use the inverse of the specified coordinate transformation to rotate
//   these unit vectors. Loop over each face, call face->Extent, and save
//   the maximum value.
//
//   G4VSolid::Inside
//
//   To decide if a point is inside, outside, or on the surface of the shape,
//   loop through all faces, and find the answer from face->Inside which gives
//   a value of "bestDistance" smaller than any other. While looping, if any
//   face->Inside returns kSurface, this value can be returned immediately.
//
//  EInside answer;
//  G4VCSGface *face = faces;
//  G4double best = kInfinity;
//  do
//  {
//    G4double distance;
//    EInside result = (*face)->Inside( p, kCarTolerance/2, distance );
//    if (result == kSurface) return kSurface;
//    if (distance < best)
//    {
//      best = distance;
//      answer = result;
//    }
//  } while( ++face < faces + numFaces );
//
//  return(answer);
//
//   G4VSolid::SurfaceNormal
//
//   Loop over all faces, call face->Normal, and return the normal to the face 
//   that is closest to the point.
//
//   G4VSolid::DistanceToIn(p)
//
//   Loop over all faces, invoking face->Distance with outgoing = false,
//   and save the answer that is smallest.
//
//   G4VSolid::DistanceToIn(p,v)
//
//   Loop over all faces, invoking face->Intersect with outgoing = false,
//   and save the answer that is smallest.
//
//   G4VSolid::DistanceToOut(p)
//
//   Loop over all faces, invoking face->Distance with outgoing = true,
//   and save the answer that is smallest.
//
//   G4VSolid::DistanceToOut(p,v)
//
//   Loop over all faces, invoking face->Intersect with outgoing = true,
//   and save the answer that is smallest. If there is more than one answer,
//   or if allBehind is false for the one answer, return validNorm as false.

// Author: David C. Williams (UCSC), 1998 - Created
// --------------------------------------------------------------------
#ifndef G4VCSGFACE_HH
#define G4VCSGFACE_HH

#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "geomdefs.hh"
#include "G4VSolid.hh"

class G4VoxelLimits;
class G4AffineTransform;
class G4SolidExtentList;

/**
 * @brief G4VCSGface is virtual base class, representing one side (or face)
 * of a CSG-like solid. It should be possible to build a CSG entirely
 * out of connecting CSG faces. Each face has an inside and outside surface,
 * the former represents the inside of the volume, the latter, the outside.
 */

class G4VCSGface
{
  public:

    /**
     * Default Constructor and Destructor.
     */
    G4VCSGface() = default;
    virtual ~G4VCSGface() = default;
  
    /**
     * Determines the distance along a line to the face.
     *  @param[in] p Position.
     *  @param[in] v Direction (assumed to be a unit vector).
     *  @param[in] outgoing Flag true, to consider only inside surfaces;
     *             false, to consider only outside surfaces.
     *  @param[in] surfTolerance Minimum distance from the surface.
     *  @param[out] distance Distance to intersection.
     *  @param[out] distFromSurface Distance from surface (along surface normal),
     *              < 0 if the point is in front of the surface.
     *  @param[out] normal Normal of surface at intersection point.
     *  @param[out] allBehind Flag, true, if entire surface is behind normal.
     *  @returns true if there is an intersection, false otherwise.
     */
    virtual G4bool Intersect( const G4ThreeVector& p, const G4ThreeVector& v,  
                              G4bool outgoing, G4double surfTolerance,
                              G4double& distance, G4double& distFromSurface,
                              G4ThreeVector& normal, G4bool& allBehind ) = 0;

    /**
     * Determines the distance of a point from either the inside or outside
     * surfaces of the face.
     *  @param[in] p Position.
     *  @param[in] outgoing Flag, true, to consider only inside surfaces
     *             or false, to consider only outside surfaces.
     *  @returns The distance to the closest surface satisfying requirements
     *           or kInfinity if no such surface exists.
     */
    virtual G4double Distance( const G4ThreeVector& p, G4bool outgoing ) = 0;
  
    /**
     * Determines whether a point is inside, outside, or on the surface of
     * the face.
     *  @param[in] p Position.
     *  @param[in] tolerance Tolerance defining the bounds of the "kSurface",
     *             nominally equal to kCarTolerance/2.
     *  @param[out] bestDistance Distance to the closest surface (in or out).
     *  @returns kInside if the point is closest to the inside surface;
     *           kOutside if the point is closest to the outside surface;
     *           kSurface if the point is withing tolerance of the surface.
     */
    virtual EInside Inside( const G4ThreeVector& p, G4double tolerance, 
                            G4double* bestDistance ) = 0;
    
    /**
     * Returns the normal of surface closest to the point.
     *  @param[in] p Position.
     *  @param[out] bestDistance Distance to the closest surface (in or out).
     *  @returns The normal of the surface nearest the point.
     */
    virtual G4ThreeVector Normal( const G4ThreeVector& p,
                                  G4double* bestDistance ) = 0;

    /**
     * Returns the face extent along the axis.
     *  @param[in] axis Unit vector defining the direction.
     *  @returns The largest point along the given axis of the face's extent.
     */
    virtual G4double Extent( const G4ThreeVector axis ) = 0;
  
    /**
     * Calculates the extent of the face for the voxel navigator.
     *  @param[in] axis The axis in which to check the shapes 3D extent against.
     *  @param[in] voxelLimit Limits along x, y, and/or z axes.
     *  @param[in] tranform A coordinate transformation on which to apply to
     *             the shape before testing.
     *  @param[out] extentList The list of (voxel) extents along the axis.
     */
    virtual void CalculateExtent( const EAxis axis, 
                                  const G4VoxelLimits& voxelLimit,
                                  const G4AffineTransform& tranform,
                                  G4SolidExtentList& extentList       ) = 0;

    /**
     * Method invoked by the copy constructor or the assignment operator.
     * Its purpose is to return a pointer to a duplicate copy of the face.
     */
    virtual G4VCSGface* Clone() = 0;

    /**
     * Returning an estimation of the face surface area, in internal units.
     */
    virtual G4double SurfaceArea() = 0;

    /**
     * Auxiliary method for GetPointOnSurface().
     */
    virtual G4ThreeVector GetPointOnFace() = 0;
};

#endif
