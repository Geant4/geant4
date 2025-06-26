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
// G4ParameterisedNavigation
//
// Class description:
//
// Utility for navigation in volumes containing a single G4PVParameterised
// volume for which voxels for the replicated volumes have been constructed.
// [Voxels MUST be along one axis only: NOT refined]

// Author: Paul Kent (CERN), August 1996
// --------------------------------------------------------------------
#ifndef G4PARAMETERISEDNAVIGATION_HH
#define G4PARAMETERISEDNAVIGATION_HH 1

#include "G4Types.hh"

#include <vector>

#include "G4VoxelNavigation.hh"
#include "G4NavigationHistory.hh"
#include "G4AffineTransform.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VPVParameterisation.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4ThreeVector.hh"
#include "G4BlockingList.hh"

/**
 * @brief G4ParameterisedNavigation is a concrete utility class for navigation
 * in volumes containing a single G4PVParameterised volume for which voxels for
 * the replicated volumes have been constructed.
 * @note Voxels MUST be along one axis only: NOT refined.
 */

class G4ParameterisedNavigation : public G4VoxelNavigation
{
  public:

    /**
     * Constructor and default Destructor.
     */
    G4ParameterisedNavigation();
    ~G4ParameterisedNavigation() override;

    /**
     * Locates voxel node based on given point. If no parameterisation axis
     * is specified, adopt default location strategy as for placements.
     *  @param[in] pHead Pointer to header of nodes to look through.
     *  @param[in] localPoint Local point
     *  @returns Pointer to the node where the given point is located.
     */
    inline G4SmartVoxelNode* ParamVoxelLocate( G4SmartVoxelHeader* pHead,
                                         const G4ThreeVector& localPoint );

    /**
     * Searches positioned volumes in mother at current top level of @p history
     * for volume containing @p globalPoint. Do not test against @p blockedVol.
     * If a containing volume is found, push it onto navigation history state.
     *  @param[in,out] history Navigation history.
     *  @param[in,out] blockedVol Blocked volume to be ignored in queries.
     *  @param[in,out] blockedNum Copy number for blocked replica volumes.
     *  @param[in,out] globalPoint Point in global coordinates system.
     *  @param[in,out] globalDirection Global direction vector.
     *  @param[in] pLocatedOnEdge Flag specifying if point is located on edge.
     *  @param[in,out] localPoint Point in local coordinates system.
     *  @returns Whether a containing volume has been found.
     */
    G4bool LevelLocate( G4NavigationHistory& history,
                  const G4VPhysicalVolume* blockedVol,
                  const G4int blockedNum,
                  const G4ThreeVector& globalPoint,
                  const G4ThreeVector* globalDirection,
                  const G4bool pLocatedOnEdge, 
                        G4ThreeVector& localPoint ) override;

    /**
     * Computes the length of a step to the next boundary.
     * Does not test against @p pBlockedPhysical. Identifies the next candidate
     * volume (if a daughter of the current volume), and returns it in:
     * pBlockedPhysical, blockedReplicaNo.
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
                                G4VPhysicalVolume *(*pBlockedPhysical),
                                G4int& blockedReplicaNo ) override;

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
    G4double ComputeSafety( const G4ThreeVector& localPoint,
                            const G4NavigationHistory& history,
                            const G4double pProposedMaxLength=DBL_MAX ) override;

    /**
     * Updates internal navigation state to take into account that location
     * has been moved, but remains within the @p motherPhysical volume.
     *  @param[in] motherPhysical Current physical volume.
     *  @param[in] localPoint Local point.
     */
    void RelocateWithinVolume( G4VPhysicalVolume* motherPhysical,
                               const G4ThreeVector& localPoint ) override;

  private:

    /**
     * Computes safety from specified point to voxel boundaries using already
     * located point.
     *  @param[in] localPoint Local point.
     *  @param[in] pAxis Axis of parameterisation.
     *  @returns Safety length from current point to voxel boundary.
     */
    G4double ComputeVoxelSafety( const G4ThreeVector& localPoint,
                                 const EAxis pAxis ) const;

    /**
     * Finds the next voxel from the current voxel and point in the specified
     * direction.
     *  @param[in] localPoint Local point.
     *  @param[in] localDirection Direction along which compute the distance.
     *  @param[in] currentStep Current step size.
     *  @param[in] pAxis Axis of parameterisation.
     *  @returns false if all voxels considered
     *           [current Step ends inside same voxel or leaves all voxels]
     *           true  otherwise
     *           [the information on the next voxel is saved].
     */
    G4bool LocateNextVoxel( const G4ThreeVector& localPoint,
                            const G4ThreeVector& localDirection,
                            const G4double currentStep,
                            const EAxis pAxis );

    /**
     * Calls virtual 'Compute' methods, and copies information if nested.
     * Method necessary to resolve cases with nested parameterisations.
     *  @param[in] num Copy number of parameterisation.
     *  @param[in] apparentPhys Potentially a PhysV or PhysT.
     *  @param[in] curParam Pointer to the parameterisation algorithm.
     *  @returns A pointer to the computed parameterised solid.
     */
    inline G4VSolid* IdentifyAndPlaceSolid( G4int num,
                                     G4VPhysicalVolume* apparentPhys, 
                                     G4VPVParameterisation* curParam );

  private:

    //  Voxel Stack information (for 1D optimisation only)
    //
    EAxis fVoxelAxis = kUndefined;
    G4int fVoxelNoSlices = 0;
    G4double fVoxelSliceWidth = 0.0; 
    std::size_t fVoxelNodeNo = 0;  
    G4SmartVoxelHeader* fVoxelHeader = nullptr;
};

#include "G4ParameterisedNavigation.icc"

#endif
