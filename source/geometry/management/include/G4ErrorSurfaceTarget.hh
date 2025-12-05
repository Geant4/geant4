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
// G4ErrorSurfaceTarget
//
// Class Description:
//
// Base class for G4ErrorTarget classes that are surfaces.

// Author: Pedro Arce (CIEMAT), September 2004
// --------------------------------------------------------------------
#ifndef G4ERRORSURFACETARGET_HH
#define G4ERRORSURFACETARGET_HH

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4ErrorTanPlaneTarget.hh"
#include "G4Plane3D.hh"

/**
 * @brief G4ErrorSurfaceTarget is a base class for G4ErrorTarget classes
 * that are surfaces.
 */

class G4ErrorSurfaceTarget : public G4ErrorTanPlaneTarget
{
  public:

    /**
     * Default Constructor and Destructor.
     */
    G4ErrorSurfaceTarget() = default;
    ~G4ErrorSurfaceTarget() override = default;

    /**
     * Computes the distance from a point to the surface in a given direction.
     *  @param[in] point The point of reference.
     *  @param[in] direc The direction vector.
     *  @returns The distance value.
     */
    G4double GetDistanceFromPoint( const G4ThreeVector& point,
                                   const G4ThreeVector& direc ) const override = 0;

    /**
     * Computes the minimal distance from a point to the surface in any
     * direction.
     *  @param[in] point The point of reference.
     *  @returns The distance value.
     */
    G4double GetDistanceFromPoint( const G4ThreeVector& point ) const override = 0;

    /**
     * Computes the plane tangent to surface at a given point.
     *  @param[in] point The point of reference.
     *  @returns The tangent plane.
     */
    G4Plane3D GetTangentPlane( const G4ThreeVector& point ) const override = 0;

    /**
     * Dumps to standard output the surface parameters.
     */
    void Dump( const G4String& msg ) const override = 0;
};

#endif
