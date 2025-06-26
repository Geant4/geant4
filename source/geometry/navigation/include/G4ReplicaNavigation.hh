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
// G4ReplicaNavigation
//
// Class description:
//
// Utility for navigation in volumes containing a single G4PVParameterised
// volume for which voxels for the replicated volumes have been constructed.
// [Voxels MUST be along one axis only: NOT refined]

// Author: Paul Kent (CERN), August 1996
// --------------------------------------------------------------------
#ifndef G4REPLICANAVIGATION_HH
#define G4REPLICANAVIGATION_HH 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "G4Types.hh"
#include "G4NavigationHistory.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4BlockingList.hh"

// Required for voxel handling
//
#include "G4SmartVoxelHeader.hh"

class G4VSolid;

struct G4ExitNormal
{
  /**
   * @brief G4ExitNormal, a bucket to hold value of Normal (3-vector), Booleans
   * for calculated and leave-behind or 'validConvex', and exiting side..
   */
   
   /** Identity of 'Side' of Replicas. Used by DistanceToOut methods. */
   enum  ESide {kNull,kRMin,kRMax,kSPhi,kEPhi,kPX,kMX,kPY,kMY,kPZ,kMZ,kMother};

   G4ThreeVector exitNormal;
   G4bool        calculated;   // Normal
   G4bool        validConvex;  // Solid locally convex
   ESide         exitSide;

 public:

   G4ExitNormal(const G4ThreeVector& norm = G4ThreeVector(0.,0.,0.),
                G4bool        calc = false,
                G4bool        valid= false,
                ESide         side = kNull )
   { exitNormal= norm; calculated= calc; validConvex=valid; exitSide=side;}
};

/**
 * @brief G4ReplicaNavigation is a utility class for navigation in volumes
 * containing a single G4PVParameterised volume for which voxels for the
 * replicated volumes have been constructed.
 * @note Voxels MUST be along one axis only: NOT refined.
 */

class G4ReplicaNavigation
{
  public:

    /**
     * Constructor and default Destructor.
     */
    G4ReplicaNavigation();
   ~G4ReplicaNavigation() = default;

    /**
     * Searches positioned volumes in mother at current top level of @p history
     * for volume containing @p globalPoint. Do not test against @p blockedVol.
     * If a containing volume is found, push it onto navigation history state.
     *  @param[in,out] history Navigation history.
     *  @param[in,out] blockedVol Blocked volume to be ignored in queries.
     *  @param[in,out] blockedNum Copy number for blocked replica volumes.
     *  @param[in,out] globalPoint Point in global coordinates system.
     *  @param[in,out] globalDirection Pointer to global direction or null.
     *  @param[in,out] localPoint Point in local coordinates system.
     *  @returns Whether a containing volume has been found.
     */
    inline G4bool LevelLocate( G4NavigationHistory& history,
                         const G4VPhysicalVolume* blockedVol,
                         const G4int blockedNum,
                         const G4ThreeVector& globalPoint,
                         const G4ThreeVector* globalDirection,
                         const G4bool pLocatedOnEdge, 
                               G4ThreeVector& localPoint );

    /**
     * Computes the length of a step to the next boundary.
     * Do not test against @p pBlockedPhysical. Identify the next candidate
     * volume (if a daughter of the current volume), and return it in:
     * pBlockedPhysical, blockedReplicaNo.
     *  @param[in] globalPoint Global point.
     *  @param[in] globalDirection Global direction vector.
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
    G4double ComputeStep( const G4ThreeVector& globalPoint,
                          const G4ThreeVector& globalDirection,
                          const G4ThreeVector& localPoint,
                          const G4ThreeVector& localDirection,
                          const G4double currentProposedStepLength,
                                G4double& newSafety,
                                G4NavigationHistory& history,
                                G4bool& validExitNormal,
                                G4bool& calculatedExitNormal,
                                G4ThreeVector &exitNormal,
                                G4bool& exiting,
                                G4bool& entering,
                                G4VPhysicalVolume* (*pBlockedPhysical),
                                G4int &blockedReplicaNo );

    /**
     * Calculates the isotropic distance to the nearest boundary from the
     * specified point in the local/global coordinate system. 
     * The localpoint utilised must be within the current volume.
     *  @param[in] globalPoint Global point.
     *  @param[in] localPoint Local point.
     *  @param[in] history Navigation history.
     *  @param[in] pProposedMaxLength Maximum step length beyond which volumes
     *             need not be checked.
     *  @returns Length from current point to closest surface.
     */
    G4double ComputeSafety( const G4ThreeVector& globalPoint,
                            const G4ThreeVector& localPoint,
                            const G4NavigationHistory& history,
                            const G4double pProposedMaxLength = DBL_MAX ) const;

    /**
     * Locates the specified point in the local/global coordinate system. 
     * The localpoint utilised must be within the current volume.
     *  @param[in] history Navigation history.
     *  @param[in] globalPoint Global point.
     *  @param[in] localPoint Local point.
     *  @param[in,out] exiting Flag to indicate whether exiting a volume.
     *  @param[in,out] notKnownInside Flag to indicate whether exiting replica.
     *  @returns If point is located inside the current replica volume.
     */
    EInside BackLocate( G4NavigationHistory &history,
                  const G4ThreeVector& globalPoint,
                        G4ThreeVector& localPoint,
                  const G4bool& exiting,
                        G4bool& notKnownInside ) const;

    /**
     * Setups transformation and transform point into local system.
     */
    void ComputeTransformation( const G4int replicaNo,
                                      G4VPhysicalVolume* pVol,
                                      G4ThreeVector& point ) const;

    /**
     * Setups transformation into local system.
     */
    void ComputeTransformation( const G4int replicaNo,
                                      G4VPhysicalVolume* pVol ) const; 

    /**
     * Computes if local point point is inside the reference volume or not.
     *  @param[in] pVol Pointer to the reference volume.
     *  @param[in] replicaNo Replica number of current volume.
     *  @param[in] localPoint Point in local coordinates system.
     *  @returns If point is inside referenced replica volume or not.
     */
    EInside Inside( const G4VPhysicalVolume* pVol,
                    const G4int replicaNo,
                    const G4ThreeVector& localPoint ) const;

    /**
     * Estimates the isotropic distance to exit the reference volume from
     * the provided local point.
     *  @param[in] pVol Pointer to the reference volume.
     *  @param[in] replicaNo Replica number of current volume.
     *  @param[in] localPoint Point in local coordinates system.
     *  @returns The isotropic distance to exit.
     */
    G4double DistanceToOut( const G4VPhysicalVolume* pVol,
                            const G4int replicaNo,
                            const G4ThreeVector& localPoint ) const;

    /**
     * Calculates the distance to exit the reference volume from the provided
     * local point, given its direction.
     *  @param[in] pVol Pointer to the reference volume.
     *  @param[in] replicaNo Replica number of current volume.
     *  @param[in] localPoint Point in local coordinates system.
     *  @param[in] localDirection The local direction.
     *  @param[in,out] candidateNormal The candidate exit notmal.
     *  @returns The distance to exit.
     */
    G4double DistanceToOut( const G4VPhysicalVolume* pVol,
                            const G4int replicaNo,
                            const G4ThreeVector& localPoint,
                            const G4ThreeVector& localDirection,
                                  G4ExitNormal& candidateNormal ) const;

    /**
     * Verbosity control.
     *  @note If level>0 && G4VERBOSE, printout can occur.
     */
    inline G4int GetVerboseLevel() const;
    inline void  SetVerboseLevel(G4int level);

    /**
     * Run navigation in "check-mode", therefore using additional
     * verifications and more strict correctness conditions.
     *  @note Is effective only with G4VERBOSE set.
     *  @param[in] mode Flag to enable/disable check-mode.
     */
    inline void CheckMode(G4bool mode);

  private:

    /**
     * Locates voxel node based on given point.
     *  @param[in] pHead Pointer to header of nodes to look through.
     *  @param[in] localPoint Local point
     *  @param[in] blocked Flag to indicate if volume is blocked or not.
     *  @returns Pointer to the node where the given point is located.
     */
    inline G4int VoxelLocate( const G4SmartVoxelHeader* pHead,
                              const G4ThreeVector& localPoint,
                              const G4int blocked=-1 ) const;

    /**
     * Computes distance to exit for phi replica.
     */
    G4double DistanceToOutPhi( const G4ThreeVector& localPoint,
                               const G4ThreeVector& localDirection,
                               const G4double width,
                               G4ExitNormal& foundNormal ) const;

    /**
     * Computes distance to exit for radial replica.
     */
    G4double DistanceToOutRad( const G4ThreeVector& localPoint,
                               const G4ThreeVector& localDirection,
                               const G4double width,
                               const G4double offset,
                               const G4int replicaNo,
                                     G4ExitNormal& foundNormal ) const;

    /**
     * Sets phi rotation of target volume.
     */
    inline void SetPhiTransformation( const G4double ang,
                                      G4VPhysicalVolume* pVol = nullptr ) const;

  private:

    // Invariants - unaltered during navigation -------------------------------
    // **********

    /** Check mode. */
    G4bool fCheck = false; 
 
    /** Verbosity. */
    G4int  fVerbose = 0;

    /** Local cached constants. */
    G4double kCarTolerance, kRadTolerance, kAngTolerance,
             halfkCarTolerance, halfkRadTolerance, halfkAngTolerance,
             fMinStep;
};

#include "G4ReplicaNavigation.icc"

#endif
