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
// G4ErrorPlaneSurfaceTarget
//
// Class Description:
//
// Limits step when track reaches a plane surface.

// Author: Pedro Arce (CIEMAT), September 2004
// --------------------------------------------------------------------
#ifndef G4ERRORPLANESURFACETARGET_HH
#define G4ERRORPLANESURFACETARGET_HH

#include "globals.hh"
#include "G4ErrorSurfaceTarget.hh"
#include "G4ThreeVector.hh"
#include "G4Normal3D.hh"
#include "G4Plane3D.hh"
#include "G4Point3D.hh"

/**
 * @brief G4ErrorPlaneSurfaceTarget is an utility class for limiting
 * the step when a track reaches a plane surface.
 */

class G4ErrorPlaneSurfaceTarget : public G4ErrorSurfaceTarget, G4Plane3D
{
  public:

    /**
     * Constructor for G4ErrorPlaneSurfaceTarget. It constructs a plane
     * by parameters: ax+by+cz+d = 0.
     */
    G4ErrorPlaneSurfaceTarget(G4double a=0., G4double b=0.,
                              G4double c=0., G4double d=0.);

    /**
     * Constructor for G4ErrorPlaneSurfaceTarget. It constructs a plane
     * by point 'p' and normal 'n'.
     */
    G4ErrorPlaneSurfaceTarget(const G4Normal3D& n,
                              const G4Point3D& p);

    /**
     * Constructor for G4ErrorPlaneSurfaceTarget. It constructs a plane
     * by three points, 'p1', 'p2', 'p3'.
     */
    G4ErrorPlaneSurfaceTarget(const G4Point3D& p1,
                              const G4Point3D& p2,
                              const G4Point3D& p3);

    /**
     * Default Destructor.
     */
    ~G4ErrorPlaneSurfaceTarget() override = default;

    /**
     * Intersects the surface with the line given by point and direction.
     *  @param[in] point The point of reference.
     *  @param[in] direc The direction vector.
     *  @returns The intersection point.
     */
    G4ThreeVector Intersect( const G4ThreeVector& point,
                             const G4ThreeVector& direc ) const;
  
    /**
     * Computes the distance from a point to the surface in a given direction.
     *  @param[in] point The point of reference.
     *  @param[in] direc The direction vector.
     *  @returns The distance value.
     */
    G4double GetDistanceFromPoint( const G4ThreeVector& point,
                                   const G4ThreeVector& direc ) const override;

    /**
     * Computes the minimal distance from a point to surface.
     *  @param[in] point The point of reference.
     *  @returns The distance value.
     */
    G4double GetDistanceFromPoint( const G4ThreeVector& pt ) const override;

    /**
     * Computes the plane tangent to itself at a given point.
     *  @param[in] point The point of reference.
     *  @returns The tangent plane.
     */
    G4Plane3D GetTangentPlane( const G4ThreeVector& point ) const override;

    /**
     * Dumps to standard output the surface parameters.
     */
    void Dump( const G4String& msg ) const override;
};

#endif
