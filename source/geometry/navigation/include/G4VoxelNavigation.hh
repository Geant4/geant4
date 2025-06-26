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
// G4VoxelNavigation
//
// Class description:
//
// Utility for navigation in volumes containing only G4PVPlacement
// daughter volumes for which voxels have been constructed.

// Author: Paul Kent (CERN), August 1996
// --------------------------------------------------------------------
#ifndef G4VOXELNAVIGATION_HH
#define G4VOXELNAVIGATION_HH 1

#include "geomdefs.hh"
#include "G4VNavigation.hh"
#include "G4NavigationHistory.hh"
#include "G4NavigationLogger.hh"
#include "G4AffineTransform.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4ThreeVector.hh"

#include "G4BlockingList.hh"

class G4VoxelSafety; 

// Required for inline implementation
//
#include "G4AuxiliaryNavServices.hh"

// Required for voxel handling & voxel stack
//
#include <vector>
#include "G4SmartVoxelProxy.hh"
#include "G4SmartVoxelNode.hh"
#include "G4SmartVoxelHeader.hh"

/**
 * @brief G4VoxelNavigation is a concrete utility class for navigation in
 * volumes containing only G4PVPlacement daughter volumes for which voxels
 * have been constructed.
 */

class G4VoxelNavigation : public G4VNavigation
{
  public:

    /**
     * Constructor and Destructor.
     */
    G4VoxelNavigation();
    ~G4VoxelNavigation() override;

    /**
     * Locates voxel node based on given point.
     *  @param[in] pHead Pointer to header of nodes to look through.
     *  @param[in] localPoint Local point
     *  @returns Pointer to the node where the given point is located.
     */
    inline G4SmartVoxelNode* VoxelLocate( G4SmartVoxelHeader* pHead,
                                    const G4ThreeVector& localPoint );

    /**
     * Searches positioned volumes in mother at current top level of @p history
     * for volume containing @p globalPoint. Do not test against @p blockedVol.
     * If a containing volume is found, push it onto navigation history state.
     *  @param[in,out] history Navigation history.
     *  @param[in,out] blockedVol Blocked volume to be ignored in queries.
     *  @param[in,out] blockedNum Copy number for blocked replica volumes.
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
                                G4VPhysicalVolume* (*pBlockedPhysical),
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
    G4double ComputeSafety( const G4ThreeVector& localpoint,
                            const G4NavigationHistory& history,
                            const G4double pMaxLength = DBL_MAX ) override;

    /**
     * Updates internal navigation state to take into account that location
     * has been moved, but remains within the @p motherPhysical volume.
     *  @param[in] motherPhysical Current physical volume.
     *  @param[in] localPoint Local point.
     */
    void RelocateWithinVolume( G4VPhysicalVolume* motherPhysical,
                               const G4ThreeVector& localPoint ) override;

    /**
     * Verbosity control.
     *  @note If level>0 && G4VERBOSE, printout can occur.
     */
    inline G4int GetVerboseLevel() const override;
    void  SetVerboseLevel(G4int level) override;

    /**
     * Enables best-possible evaluation of isotropic safety.
     */
    inline void EnableBestSafety( G4bool flag = false );

  protected:

    /**
     * Computes safety from specified point to voxel boundaries using already
     * located point.
     *  @param[in] localPoint Local point.
     *  @returns Safety length from current point to voxel boundary.
     */
    G4double ComputeVoxelSafety( const G4ThreeVector& localPoint ) const;

    /**
     * Finds the next voxel from the current voxel and point in the specified
     * direction.
     *  @param[in] localPoint Local point.
     *  @param[in] localDirection Direction along which compute the distance.
     *  @param[in] currentStep Current step size.
     *  @returns false if all voxels considered
     *           [current Step ends inside same voxel or leaves all voxels]
     *           true  otherwise
     *           [the information on the next voxel is saved].
     */
    G4bool LocateNextVoxel( const G4ThreeVector& localPoint,
                            const G4ThreeVector& localDirection,
                            const G4double currentStep );

  private:  // Logging functions

    void PreComputeStepLog  (const G4VPhysicalVolume* motherPhysical,
                                   G4double motherSafety,
                             const G4ThreeVector& localPoint);
    void AlongComputeStepLog(const G4VSolid* sampleSolid,
                             const G4ThreeVector& samplePoint,
                             const G4ThreeVector& sampleDirection,
                             const G4ThreeVector& localDirection,
                                   G4double sampleSafety,
                                   G4double sampleStep);
    void PostComputeStepLog (const G4VSolid* motherSolid,
                             const G4ThreeVector& localPoint,
                             const G4ThreeVector& localDirection,
                                   G4double motherStep,
                                   G4double motherSafety);
    void ComputeSafetyLog   (const G4VSolid* solid,
                             const G4ThreeVector& point,
                                   G4double safety,
                                   G4bool banner);
    inline void PrintDaughterLog (const G4VSolid* sampleSolid,
                                  const G4ThreeVector& samplePoint,
                                        G4double sampleSafety,
                                        G4double sampleStep);   
  protected:

    /** Blocked volumes. */
    G4BlockingList fBList;

    // -----------------------------------------------------------------------
    // BEGIN Voxel Stack information

    /** Voxels depth.
     *   @note fVoxelDepth==0+ => fVoxelAxisStack(0+) contains axes of voxel
     *         fVoxelDepth==-1 -> not in voxel.
     */
    G4int fVoxelDepth = -1;

    /** Voxel axes. */
    std::vector<EAxis> fVoxelAxisStack;

    /** No slices per voxel at each level. */
    std::vector<G4int> fVoxelNoSlicesStack;

    /** Width of voxels at each level. */
    std::vector<G4double> fVoxelSliceWidthStack; 

    /** Node no point is inside at each level. */
    std::vector<G4int> fVoxelNodeNoStack;    

    /** Voxel headers at each level. */
    std::vector<G4SmartVoxelHeader*> fVoxelHeaderStack;

    /** Node containing last located point. */
    G4SmartVoxelNode* fVoxelNode = nullptr;

    // END Voxel Stack information
    // -----------------------------------------------------------------------

    /** Helper object for Voxel Safety. */
    G4VoxelSafety* fpVoxelSafety = nullptr;

    /** Surface tolerance. */
    G4double fHalfTolerance;

    /** Flag for best safety. */
    G4bool fBestSafety = false; 

    /** Verbosity logger. */
    G4NavigationLogger* fLogger;
};

#include "G4VoxelNavigation.icc"

#endif
