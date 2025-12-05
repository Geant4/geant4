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
// G4ErrorCylSurfaceTarget
//
// Class Description:
//
// Limits step when track reaches a cylindrical surface.

// Author: Pedro Arce (CIEMAT), September 2004
// --------------------------------------------------------------------
#ifndef G4ERRORCYLSURFACETARGET_HH
#define G4ERRORCYLSURFACETARGET_HH

#include "globals.hh"
#include "G4ErrorSurfaceTarget.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4Plane3D.hh"

/**
 * @brief G4ErrorCylSurfaceTarget is a utility class for limiting
 * the step when a track reaches a cylindrical surface.
 */

class G4ErrorCylSurfaceTarget : public G4ErrorSurfaceTarget
{
  public:

    /**
     * Constructor for G4ErrorCylSurfaceTarget. It constructs a cylindrical
     * surface by radius, translation and rotation.
     *  @param[in] radius Cylinder radius.
     *  @param[in] trans Translation vector.
     *  @param[in] rotm Rotation matrix.
     */
    G4ErrorCylSurfaceTarget( const G4double& radius,
                             const G4ThreeVector& trans = G4ThreeVector(),
                             const G4RotationMatrix& rotm = G4RotationMatrix() );

    /**
     * Constructor for G4ErrorCylSurfaceTarget. It constructs a cylindrical
     * surface by radius and affine transformation.
     *  @param[in] radius Cylinder radius.
     *  @param[in] trans The affine transformation in input.
     */
    G4ErrorCylSurfaceTarget( const G4double& radius,
                             const G4AffineTransform& trans );

    /**
     * Default Destructor.
     */
    ~G4ErrorCylSurfaceTarget() override = default;

    /**
     * Intersects the cylindrical surface with the line in local (cylinder)
     * coordinates, given by point and direction.
     *  @param[in] point The point of reference.
     *  @param[in] direc The direction vector.
     *  @returns The intersection point.
     */
    G4ThreeVector IntersectLocal( const G4ThreeVector& point,
                                  const G4ThreeVector& direc ) const;

    /**
     * Computes the distance from a point to the cylindrical surface in a
     * given direction.
     *  @param[in] point The point of reference.
     *  @param[in] direc The direction vector.
     *  @returns The distance value.
     */
    G4double GetDistanceFromPoint( const G4ThreeVector& point,
                                   const G4ThreeVector& direc ) const override;

    /**
     * Computes the minimal distance from a point to the cylindrical surface
     * in any direction.
     *  @param[in] point The point of reference.
     *  @returns The distance value.
     */
    G4double GetDistanceFromPoint( const G4ThreeVector& point ) const override;

    /**
     * Computes the plane tangent to cylindrical surface at a given point.
     *  @param[in] point The point of reference.
     *  @returns The tangent plane.
     */
    G4Plane3D GetTangentPlane( const G4ThreeVector& point ) const override;

    /**
     * Dumps to standard output the cylindrical surface parameters.
     */
    void Dump( const G4String& msg ) const override;

  private:

     G4double fradius;
     G4AffineTransform ftransform;
};

#endif
