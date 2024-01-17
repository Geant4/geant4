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
// class G4VNavigation
//
// Class description:
//
// Navigation interface common between all navigator types.

// Author: G. Amadio - CERN, March 2022
// --------------------------------------------------------------------
#ifndef G4VNAVIGATION_HH
#define G4VNAVIGATION_HH

#include "G4ThreeVector.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4NavigationHistory;

/**
 * @brief G4VNavigation class holds the common navigation interface
 * for all geometry navigator types.
 */

class G4VNavigation
{
 public:
  /** Virtual Destructor. */
  virtual ~G4VNavigation() {}

  /**
   * Search positioned volumes in mother at current top level of @p history
   * for volume containing @p globalPoint. Do not test against @p blockedVol.
   * If a containing volume is found, push it onto navigation history state.
   * @param[in,out] history Navigation history.
   * @param[in,out] blockedVol Blocked volume that should be ignored in queries.
   * @param[in,out] blockedNum Copy number for blocked replica volumes.
   * @param[in,out] globalPoint Global point
   * @param[in,out] globalDirection Pointer to global direction or null pointer.
   * @param[in,out] localPoint = global point in local system on entry, point
   *                in new system on exit.
   * @returns Whether a containing volume has been found.
   */
  virtual G4bool LevelLocate(G4NavigationHistory& history,
                             const G4VPhysicalVolume* blockedVol,
                             const G4int blockedNum,
                             const G4ThreeVector& globalPoint,
                             const G4ThreeVector* globalDirection,
                             const G4bool pLocatedOnEdge,
                             G4ThreeVector& localPoint) = 0;

  /**
   * Compute the length of a step to the next boundary.
   * Do not test against @p pBlockedPhysical. Identify the next candidate volume
   * (if a daughter of current volume), and return it in pBlockedPhysical,
   * blockedReplicaNo.
   * @param[in] localPoint Local point
   * @param[in] localDirection Pointer to local direction or null pointer.
   * @param[in] currentProposedStepLength Current proposed step length.
   * @param[in,out] newSafety New safety.
   * @param[in,out] history Navigation history.
   * @param[in,out] validExitNormal Flag to indicate whether exit normal is
   * valid or not.
   * @param[in,out] exitNormal Exit normal.
   * @param[in,out] entering Flag to indicate whether we are entering a volume.
   * @param[in,out] exiting Flag to indicate whether we are exiting a volume.
   * @param[in,out] pBlockedPhysical Blocked physical volume that should be
   * ignored in queries.
   * @param[in,out] blockedReplicaNo Copy number for blocked replica volumes.
   * @returns Length from current point to next boundary surface along @p
   * localDirection.
   */
  virtual G4double ComputeStep(const G4ThreeVector& localPoint,
                               const G4ThreeVector& localDirection,
                               const G4double currentProposedStepLength,
                               G4double& newSafety,
                               G4NavigationHistory& history,
                               G4bool& validExitNormal,
                               G4ThreeVector& exitNormal,
                               G4bool& exiting,
                               G4bool& entering,
                               G4VPhysicalVolume*(*pBlockedPhysical),
                               G4int& blockedReplicaNo) = 0;

  /**
   * Compute the distance to the closest surface.
   * @param[in] globalPoint Global point.
   * @param[in] history Navigation history.
   * @param[in] pMaxLength Maximum step length beyond which volumes need not be
   * checked.
   * @returns Length from current point to closest surface.
   */
  virtual G4double ComputeSafety(const G4ThreeVector& globalpoint,
                                 const G4NavigationHistory& history,
                                 const G4double pMaxLength = DBL_MAX) = 0;

  /**
   * Update internal navigation state to take into account that location
   * has been moved, but remains within the @p motherPhysical volume.
   *  @param[in] motherPhysical Current physical volume.
   *  @param[in] localPoint Local point.
   */
  virtual void RelocateWithinVolume(G4VPhysicalVolume* /* motherPhysical */,
                                    const G4ThreeVector& /* localPoint */)
  {
    /* do nothing by default */
  }

  /** Get current verbosity level */
  virtual G4int GetVerboseLevel() const { return fVerbose; }

  /** Set current verbosity level */
  virtual void SetVerboseLevel(G4int level) { fVerbose = level; }

  /**
   * Set check mode.
   * When enabled, forces navigator to run in "check mode", hence using
   * additional verifications and stricter condictions for ensuring correctness.
   * Effective only when G4VERBOSE is enabled.
   */
  void CheckMode(G4bool mode) { fCheck = mode; }

 protected:
  G4int fVerbose = 0;
  G4bool fCheck = false;
};

#endif
