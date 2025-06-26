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
// G4NormalNavigation
//
// Class description:
//
// Utility for navigation in volumes containing only G4PVPlacement
// daughter volumes.

// Author: Paul Kent (CERN), August 1996
// --------------------------------------------------------------------
#ifndef G4NORMALNAVIGATION_HH
#define G4NORMALNAVIGATION_HH 1

#include <iomanip>

#include "G4VNavigation.hh"
#include "G4NavigationHistory.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4ThreeVector.hh"
#include "G4AuxiliaryNavServices.hh"

class G4NavigationLogger;

/**
 * @brief G4NormalNavigation is a concrete utility class for navigation in
 * volumes containing only G4PVPlacement daughter volumes.
 */

class G4NormalNavigation : public G4VNavigation
{
  public:

    /**
     * Constructor and Destructor.
     */
    G4NormalNavigation();
    ~G4NormalNavigation() override;

    /**
     * Searches positioned volumes in mother at current top level of history
     * for volume containing @p globalPoint. Do not test against @p blockedVol.
     * If a containing volume is found, push it onto navigation history state
     * and return true, else return false (the point lying in the mother but
     * not any of the daughters).
     *  @param[in,out] history Navigation history.
     *  @param[in,out] blockedVol Blocked volume to be ignored in queries.
     *  @param[in,out] blockedNum Copy number for blocked volumes.
     *  @param[in,out] globalPoint Point in global coordinates system.
     *  @param[in,out] globalDirection Pointer to global direction or null.
     *  @param[in] pLocatedOnEdge Flag specifying if point is located on edge.
     *  @param[in,out] localPoint Point in local coordinates system.
     *  @returns Whether a containing volume has been found.
     */
    inline G4bool LevelLocate( G4NavigationHistory& history,
                  const G4VPhysicalVolume* blockedVol,
                  const G4int blockedNum,
                  const G4ThreeVector& globalPoint,
                  const G4ThreeVector* globalDirection,
                  const G4bool pLocatedOnEdge, 
                        G4ThreeVector& localPoint) final;

    /**
     * Computes the length of a step to the next boundary.
     * Does not test against @p pBlockedPhysical. Identifies the next candidate
     * volume (if a daughter of the current volume), and returns it in:
     * pBlockedPhysical, blockedReplicaNo.
     *  @param[in] localPoint Local point.
     *  @param[in] localDirection Pointer to local direction or null pointer.
     *  @param[in] currentProposedStepLength Current proposed step length.
     *  @param[in,out] newSafety New safety.
     *  @param[in,out] history Navigation history.
     *  @param[in,out] validExitNormal Flag to indicate whether exit normal is
     *                 valid or not.
     *  @param[in,out] exitNormal Exit normal.
     *  @param[in,out] exiting Flag to indicate whether exiting a volume.
     *  @param[in,out] entering Flag to indicate whether entering a volume.
     *  @param[in,out] pBlockedPhysical Blocked physical volume that should be
     *                 ignored in queries.
     *  @param[in,out] blockedReplicaNo Copy number for blocked replica volumes.
     *  @returns Length from current point to next boundary surface along
     *           @p localDirection.
     */
    G4double ComputeStep( const G4ThreeVector& localPoint,
                          const G4ThreeVector& localDirection,
                          const G4double currentProposedStepLength,
                                G4double& newSafety,
                                G4NavigationHistory& history,
                                G4bool& validExitNormal,
                                G4ThreeVector& exitNormal,
                                G4bool& exiting,
                                G4bool& entering,
                                G4VPhysicalVolume* (*pBlockedPhysical),
                                G4int& blockedReplicaNo ) final;

    /**
     * Calculates the isotropic distance to the nearest boundary from the
     * specified point in the local coordinate system. 
     * The localpoint utilised must be within the current volume.
     *  @param[in] localPoint Local point.
     *  @param[in] history Navigation history.
     *  @param[in] pMaxLength Maximum step length beyond which volumes
     *             need not be checked.
     *  @returns Length from current point to closest surface.
     */
    G4double ComputeSafety( const G4ThreeVector& localpoint,
                            const G4NavigationHistory& history,
                            const G4double pMaxLength=DBL_MAX ) final;

    /**
     * Verbosity control.
     *  @note If level>0 && G4VERBOSE, printout can occur.
     */
    G4int GetVerboseLevel() const final;
    void  SetVerboseLevel(G4int level) final;

  private:

    G4NavigationLogger* fLogger;
};

#include "G4NormalNavigation.icc"

#endif
