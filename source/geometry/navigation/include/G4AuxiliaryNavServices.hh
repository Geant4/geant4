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
// G4AuxiliaryNavServices
//
// Class description:
//
// Utility class for navigation.

// Author: Paul Kent (CERN), August 1996
// --------------------------------------------------------------------
#ifndef G4AuxiliaryNavServices_hh
#define G4AuxiliaryNavServices_hh 1

#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4VSolid.hh"
#include "G4AffineTransform.hh"

/**
 * @brief G4AuxiliaryNavServices a utility class for navigation.
 */

class G4AuxiliaryNavServices
{
  public:

    /**
     * Is the track (point, direction) inside the solid 'sampleSolid' ? 
     *  @param[in] sampleSolid Pointer to the shape to check.
     *  @param[in,out] localPoint Point in local coordinates system.
     *  @param[in,out] globalDirection Pointer to global direction or null.
     *  @param[in] sampleTransform Affine transformation in space.
     *  @param[in] pLocatedOnEdge Flag specifying if point is located on edge.
     * @returns True if we are going to enter the volume, which is the case
     *          if the point is inside, or the point is on the surface and
     *          the direction points inside or along it. Else returns false.
     */
    static G4bool CheckPointOnSurface( const G4VSolid* sampleSolid, 
                                       const G4ThreeVector& localPoint, 
                                       const G4ThreeVector* globalDirection, 
                                       const G4AffineTransform& sampleTransform,
                                       const G4bool locatedOnEdge);

    /**
     * Is the track (point, direction) exiting the solid 'sampleSolid' ? 
     *  @returns True if we are going to exit the volume.
     *  @param[in] sampleSolid Pointer to the shape to check.
     *  @param[in,out] localPoint Point in local coordinates system.
     *  @param[in,out] globalDirection Pointer to global direction or null.
     *  @param[in] sampleTransform Affine transformation in space.
     */
    static G4bool CheckPointExiting( const G4VSolid* sampleSolid, 
                                     const G4ThreeVector& localPoint, 
                                     const G4ThreeVector* globalDirection, 
                                     const G4AffineTransform& sampleTransform );

    /**
     * Prints global values of Cartesian, Radial and Angle Tolerances.
     */
    static void ReportTolerances();
};

#include "G4AuxiliaryNavServices.icc"

#endif
