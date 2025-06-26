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
// G4RegularNavigation
//
// Class description:
//
// Utility for fast navigation in volumes containing a regular
// parameterisation. If two contiguous voxels have the same material,
// navigation does not stop at the surface.

// Author: Pedro Arce (CIEMAT), May 2007
// --------------------------------------------------------------------
#ifndef G4RegularNavigation_HH
#define G4RegularNavigation_HH 1

#include <vector>

#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "G4VNavigation.hh"

class G4NormalNavigation;
class G4VPhysicalVolume;
class G4Navigator;
class G4NavigationHistory;

/**
 * @brief G4RegularNavigation is a concrete utility class for fast navigation
 * in volumes containing a regular parameterisation. If two contiguous voxels
 * have the same material, navigation does not stop at the surface.
 */

class G4RegularNavigation : public G4VNavigation
{
  public:
  
    /**
     * Constructor and Destructor.
     */
    G4RegularNavigation();
   ~G4RegularNavigation() override;
  
    /**
     * Locates a point using its position with respect to regular
     * parameterisation container volume.
     *  @param[in,out] history Navigation history.
     *  @param[in,out] blockedVol Blocked volume to be ignored in queries.
     *  @param[in,out] blockedNum Copy number for blocked replica volumes.
     *  @param[in,out] globalPoint Point in global coordinates system.
     *  @param[in,out] globalDirection Global direction vector.
     *  @param[in,out] localPoint Point in local coordinates system.
     *  @returns Whether a containing volume has been found.
     */
    G4bool LevelLocate(      G4NavigationHistory& history,
                       const G4VPhysicalVolume* blockedVol,
                       const G4int blockedNum,
                       const G4ThreeVector& globalPoint,
                       const G4ThreeVector* globalDirection,
                       const G4bool pLocatedOnEdge, 
                             G4ThreeVector& localPoint ) final;

    /**
     * Method never called because to be called the daughter has to be a
     * 'regular' volume. This would only happen if the track is in the
     * mother of voxels volume. But the voxels fill completely their mother,
     * so when a track enters the mother it automatically enters a voxel.
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
                                G4VPhysicalVolume *(*pBlockedPhysical),
                                G4int& blockedReplicaNo ) final;
  
    /**
     * Computes the step skipping surfaces when they separate voxels with
     * equal materials. Loops to voxels until a different material is found:
     * invokes G4NormalNavigation::ComputeStep() in each voxel and moves the
     * point to the next voxel.
     *  @param[in] localPoint Local point.
     *  @param[in] localDirection Local direction vector.
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
     *  @param[in] pCurrentPhysical Pointer to current volume.
     *  @returns Length from current point to next boundary surface along
     *           the direction.
     */
    G4double ComputeStepSkippingEqualMaterials( 
                                G4ThreeVector& localPoint,
                          const G4ThreeVector& localDirection,
                          const G4double currentProposedStepLength,
                                G4double& newSafety,
                                G4NavigationHistory& history,
                                G4bool& validExitNormal,
                                G4ThreeVector& exitNormal,
                                G4bool& exiting,
                                G4bool& entering,
                                G4VPhysicalVolume *(*pBlockedPhysical),
                                G4int& blockedReplicaNo,
                                G4VPhysicalVolume* pCurrentPhysical);

    /**
     * Method never called because to be called the daughter has to be a
     * 'regular' volume. This would only happen if the track is in the
     * mother of voxels volume. But the voxels fill completely their mother,
     * so when a track enters the mother it automatically enters a voxel.
     */
    G4double ComputeSafety( const G4ThreeVector& localPoint,
                            const G4NavigationHistory& history,
                            const G4double pProposedMaxLength = DBL_MAX ) final;

    /**
     * Setter for normal navigation.
     */
    void SetNormalNavigation( G4NormalNavigation* fnormnav );

  private:

    /** Cached pointer to normal navigation. */
    G4NormalNavigation* fnormalNav = nullptr;

    /** Surface tolerance. */
    G4double kCarTolerance;

    /** Cached minimum step. */
    G4double fMinStep;
 
    /** Whether the last ComputeStep moved Zero. Used to check for edges. */
    G4bool fLastStepWasZero = false;

    /** Number of preceding moves that were 0. Reset to 0 after finite step. */
    G4int fNumberZeroSteps = 0;

    /** After this many failed/zero steps, act (push etc). */
    G4int fActionThreshold_NoZeroSteps = 2;  

    /** After this many failed/zero steps, abandon track. */
    G4int fAbandonThreshold_NoZeroSteps = 25; 

    /** Maximum number of steps a track can travel skipping voxels
        (if there are more, track is assumed to be stuck and it is killed). */
    G4int fNoStepsAllowed = 10000;
};

#endif
