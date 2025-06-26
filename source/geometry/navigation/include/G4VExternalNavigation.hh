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
// G4VExternalNavigation
//
// Class description:
//
// Pure virtual class to be specialised by the user for tracking with
// an external navigation

// Authors: V.Vlachoudis, G.Cosmo (CERN), 2019
// --------------------------------------------------------------------
#ifndef G4VEXTERNALNAVIGATION_HH
#define G4VEXTERNALNAVIGATION_HH 1

#include "G4LogicalVolume.hh"
#include "G4NavigationHistory.hh"
#include "G4ThreeVector.hh"
#include "G4VNavigation.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"

/**
 * @brief G4VExternalNavigation is a pure virtual class to be specialised
 * by the user for tracking with an external navigation.
 */

class G4VExternalNavigation : public G4VNavigation
{
  public:

    /**
     * Constructor and Destructor.
     */
    G4VExternalNavigation();
    ~G4VExternalNavigation() override;

    /**
     * Cloning method, pure virtual.
     */
    virtual G4VExternalNavigation* Clone() = 0;

    // Optional methods - may be necessary under particular circumstances

    /**
     * Special 'Inside' call that includes direction of next motion.
     * Provided for potential optimisations.
     *  @param[in] solid Solid to be considered.
     *  @param[in] position Point to be checked.
     *  @param[in] direction Not used.
     *  @returns Whether the point is inside the solid or not.
     */
    virtual EInside Inside( const G4VSolid*      solid,
                            const G4ThreeVector& position,
                            const G4ThreeVector& direction );

    /**
     * Updates any relevant internal state to take account that the location
     * has been moved to 'localPoint' and that it remains in the current
     * (mother) physical volume 'motherPhysical'.
     *  @note Default action is do-nothing; only implemented by the concrete
     *        navigator class.
     *  @param[in] motherPhysical Volume to be considered.
     *  @param[in] localPoint Point to be checked.
     */
     void RelocateWithinVolume( G4VPhysicalVolume* motherPhysical,
                                const G4ThreeVector& localPoint ) override;
};

#endif
