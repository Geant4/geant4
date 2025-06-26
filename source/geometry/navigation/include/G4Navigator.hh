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
// G4Navigator
//
// Class description:
//
// A class for use by the tracking management, able to obtain/calculate
// dynamic tracking time information such as the distance to the next volume,
// or to find the physical volume containing a given point in the world
// reference system. The navigator maintains a transformation history and
// other information to optimise the tracking time performance.

// Original author: Paul Kent (CERN), July 1995-1996
//
// - Made Navigator Abstract                   G. Cosmo,      Nov  2003
// - Added check mode                          G. Cosmo,      Mar  2004
// - Zero step protections                     J.A. / G.C.,   Nov  2004
// --------------------------------------------------------------------
#ifndef G4NAVIGATOR_HH
#define G4NAVIGATOR_HH 1

#include "geomdefs.hh"

#include "G4ThreeVector.hh"
#include "G4AffineTransform.hh"
#include "G4RotationMatrix.hh"

#include "G4LogicalVolume.hh"             // Used in inline methods
#include "G4TouchableHandle.hh"           //    "         "

#include "G4NavigationHistory.hh"
#include "G4NormalNavigation.hh"
#include "G4VoxelNavigation.hh"
#include "G4ParameterisedNavigation.hh"
#include "G4ReplicaNavigation.hh"
#include "G4RegularNavigation.hh"
#include "G4VExternalNavigation.hh"

#include <iostream>

class G4VPhysicalVolume;
class G4SafetyCalculator;

/**
 * @brief G4Navigator is a class for use by the tracking management, able to
 * obtain/calculate dynamic tracking time information such as the distance to
 * the next volume, or to find the physical volume containing a given point in
 * the world reference system. The navigator maintains a transformation history
 * and other information to optimise the tracking time performance.
 */

class G4Navigator
{
  public:

    friend std::ostream& operator << (std::ostream &os, const G4Navigator &n);

    /**
     * Constructor - initialisers and setup.
     */
    G4Navigator();

    /**
     * Copy constructor & assignment operator not allowed.
     */
    G4Navigator(const G4Navigator&) = delete;
    G4Navigator& operator=(const G4Navigator&) = delete;

    /**
     * Destructor.
     */
    virtual ~G4Navigator();

    /**
     * Calculates the distance to the next boundary intersected along the
     * specified NORMALISED vector direction and from the specified point in
     * the global coordinate system.
     * LocateGlobalPointAndSetup() or LocateGlobalPointWithinVolume() must
     * have been called with the same global point prior to this call.
     * The isotropic distance to the nearest boundary is also calculated
     * (usually an underestimate). The current proposed Step length is used
     * to avoid intersection calculations: if it can be determined that the
     * nearest boundary is >pCurrentProposedStepLength away, kInfinity
     * is returned together with the computed isotropic safety distance.
     *  @note Geometry must be closed.
     *  @param[in] pGlobalPoint The point in global coordinates system.
     *  @param[in] pDirection The normalised vector direction.
     *  @param[in] pCurrentProposedStepLength Current proposed step length.
     *  @param[in,out] newSafety New safety.
     *  @returns Length from current point to next boundary surface along
     *           @p pDirection.
     */
    virtual G4double ComputeStep(const G4ThreeVector& pGlobalPoint,
                                 const G4ThreeVector& pDirection,
                                 const G4double pCurrentProposedStepLength,
                                       G4double& pNewSafety);

    /**
     * Same as ComputeStep() above, but does not affect/modify the state
     * of the Navigator.
     */
    G4double CheckNextStep(const G4ThreeVector& pGlobalPoint,
                           const G4ThreeVector& pDirection,
                           const G4double pCurrentProposedStepLength,
                                 G4double& pNewSafety); 

    /**
     * Resets the geometrical hierarchy and searches for the volumes deepest
     * in the hierarchy containing the point in the global coordinates space.
     * The direction is used to check if a volume is entered.
     * The search begin is the geometrical hierarchy at the location of the
     * last located point, or the endpoint of the previous Step if
     * SetGeometricallyLimitedStep() has been called immediately before.
     *  @note: In order to call this the geometry MUST be closed.
     *  @param[in] point The point in global coordinates system.
     *  @param[in] direction The normalised vector direction.
     *  @param[in] h The touchable history to be used for initialisation.
     *  @returns The pointer to the physical volume where point is located.
     */
    virtual
    G4VPhysicalVolume* ResetHierarchyAndLocate(const G4ThreeVector& point,
                                               const G4ThreeVector& direction,
                                               const G4TouchableHistory& h);

    /**
     * Searches the geometrical hierarchy for the volumes deepest in hierarchy
     * containing the point in the global coordinate space. Two main cases
     * are:
     *  i) If pRelativeSearch=false it makes use of no previous/state
     *     information. Returns the physical volume containing the point, 
     *     with all previous mothers correctly set up.
     * ii) If pRelativeSearch is set to true, the search begin is the
     *     geometrical hierarchy at the location of the last located point,
     *     or the endpoint of previous Step if SetGeometricallyLimitedStep()
     *     has been called immediately before.
     * The direction is used (to check if a volume is entered) if either
     *   - the argument ignoreDirection is false, or
     *   - the Navigator has determined that it is on an edge shared by two
     *     or more volumes (this is state information).
     *  @note In order to call this the geometry MUST be closed.
     *  @param[in] point The point in global coordinates system.
     *  @param[in] direction The normalised vector direction.
     *  @param[in] pRelativeSearch Flag to specify where search starts from.
     *  @param[in] ignoreDirection Flag to specify if to use direction or not.
     *  @returns The pointer to the physical volume where point is located.
     */
    virtual
    G4VPhysicalVolume* LocateGlobalPointAndSetup(const G4ThreeVector& point,
                                       const G4ThreeVector* direction = nullptr,
                                       const G4bool pRelativeSearch = true,
                                       const G4bool ignoreDirection = true);

    /**
     * Notifies the Navigator that a track has moved to the new Global point
     * 'position', that is known to be within the current safety.
     * No check is performed to ensure that it is within  the volume. 
     * This method can be called instead of LocateGlobalPointAndSetup() ONLY
     * if the caller is certain that the new global point (position) is inside
     * the same volume as the previous position. Usually this can be guaranteed
     * only if the point is within safety.
     *  @param[in] position The position point in global coordinates system.
     */
    virtual
    void LocateGlobalPointWithinVolume(const G4ThreeVector& position);

    /**
     * It first searches the geometrical hierarchy like the above method
     * LocateGlobalPointAndSetup(), then it uses the volume found and its
     * navigation history to update the touchable handle.
     *  @param[in] position The point in global coordinates system.
     *  @param[in] direction The normalised vector direction.
     *  @param[in,out] oldTouchableToUpdate Touchable handle to update.
     *  @param[in] RelativeSearch Flag to specify where search starts from.
     */
    inline void LocateGlobalPointAndUpdateTouchableHandle(
                  const G4ThreeVector&       position,
                  const G4ThreeVector&       direction,
                        G4TouchableHandle&   oldTouchableToUpdate,
                  const G4bool               RelativeSearch = true);
    /**
     * Same as the method above LocateGlobalPointAndUpdateTouchableHandle(),
     * except that a pointer to G4VTouchable is used for updating the touchable.
     */
    inline void LocateGlobalPointAndUpdateTouchable(
                  const G4ThreeVector&       position,
                  const G4ThreeVector&       direction,
                        G4VTouchable*        touchableToUpdate,
                  const G4bool               RelativeSearch = true);

    /**
     * Same as the method above LocateGlobalPointAndUpdateTouchable(),
     * except that direction is not specified.
     */
    inline void LocateGlobalPointAndUpdateTouchable(
                  const G4ThreeVector&       position,
                        G4VTouchable*        touchableToUpdate,
                  const G4bool               RelativeSearch = true);

    /**
     * Informs the navigator that the previous Step calculated by the
     * geometry was taken in its entirety.
     */
    inline void SetGeometricallyLimitedStep();

    /**
     * Calculates the isotropic distance to the nearest boundary from the
     * specified point in the global coordinate system. 
     *  @note The geometry must be closed.
     *  @param[in] globalpoint The point in global coordinates system.
     *             The point must be within the current volume.
     *  @param[in] pProposedMaxLength The proposed maximum length is used
     *             to avoid volume safety calculations.
     *  @param[in] keepState Flag to instruct keeping the state (default true)
     *             to ensure minimum side effects from the call.
     *  @returns Length from current point to closest boundary surface.
     *           The value returned is usually an underestimate.  
     */
    virtual G4double ComputeSafety(const G4ThreeVector& globalpoint,
                                   const G4double pProposedMaxLength = DBL_MAX,
                                   const G4bool keepState = true);

    /**
     * Returns the current  world (topmost) volume.
     */
    inline G4VPhysicalVolume* GetWorldVolume() const;

    /**
     * Sets the world (topmost) volume. This must be positioned at the
     * origin (0,0,0) and unrotated.
     */
    inline void SetWorldVolume(G4VPhysicalVolume* pWorld);

    /**
     * Touchable creation method.
     *  @note Caller has deletion responsibility.
     *  @returns A pointer to the allocated touchable history.
     */
    inline G4TouchableHistory* CreateTouchableHistory() const;

    /**
     * Touchable creation method, given a history.
     *  @note Caller has deletion responsibility.
     *  @param[in] h Pointer to a navigation history to copy from.
     *  @returns A pointer to the allocated touchable history.
     */
    inline G4TouchableHistory* CreateTouchableHistory(const G4NavigationHistory* h) const;

    /**
     * Returns a reference counted handle to a touchable history.
     */
    virtual G4TouchableHandle CreateTouchableHistoryHandle() const;

    /**
     * Obtains the Normal vector to a surface (in local coordinates)
     * pointing out of previous volume and into current volume
     * Convention: the *local* normal is in the coordinate system of the
     * *final* volume. The method takes full care about how to calculate
     * this normal, but if the surfaces are not convex it will return
     * valid=false.
     *  @note Can only be called if the Navigator's last Step has crossed a
     *        volume geometrical boundary.
     *  @note Normals are not available for replica volumes (i.e. valid=false).
     *  @param[in,out] valid Flag indicating if normal is valid.
     *  @returns A Exit Surface Normal vector and validity too.
     */
    virtual G4ThreeVector GetLocalExitNormal(G4bool* valid);

    /**
     * Obtains the Normal vector to a surface (in local coordinates)
     * pointing out of previous volume and into current volume, and
     * checks the current point against expected 'local' value.
     * Convention: the *local* normal is in the coordinate system of the
     * *final* volume. The method takes full care about how to calculate
     * this normal, but if the surfaces are not convex it will return
     * valid=false.
     *  @note Can only be called if the Navigator's last Step has crossed a
     *        volume geometrical boundary.
     *  @note Normals are not available for replica volumes (i.e. valid=false).
     *  @param[in] point Point in global coordinates system to compare to.
     *  @param[in,out] valid Flag indicating if normal is valid.
     *  @returns A Exit Surface Normal vector and validity too.
     */
    virtual G4ThreeVector GetLocalExitNormalAndCheck(const G4ThreeVector& point,
                                                           G4bool* valid);
    /**
     * Obtains the Normal vector to a surface (in global coordinates)
     * pointing out of previous volume and into current volume
     * The method takes full care about how to calculate the normal,
     * but if the surfaces are not convex it will return valid=false.
     *  @note Can only be called if the Navigator's last Step has crossed a
     *        volume geometrical boundary.
     *  @note Normals are not available for replica volumes (i.e. valid=false).
     *  @param[in] point Point in global coordinates system to compare to.
     *  @param[in,out] valid Flag indicating if normal is valid.
     *  @returns A Exit Surface Normal vector and validity too.
     */
    virtual G4ThreeVector GetGlobalExitNormal(const G4ThreeVector& point,
                                                    G4bool* valid);

    /**
     * Verbosity control.
     *  @note If level>0 && G4VERBOSE, printout can occur.
     */
    inline G4int GetVerboseLevel() const;
    inline void  SetVerboseLevel(G4int level);

    /**
     * Verify if the navigator is active.
     */
    inline G4bool IsActive() const;

    /**
     * Activate/inactivate the navigator.
     */
    inline void  Activate(G4bool flag);

    /**
     * The purpose of this function is to inform the caller if the track is
     * entering a daughter volume while exiting from the current volume.
     *  @note It is not guaranteed to work if SetGeometricallyLimitedStep()
     *        was not called when it should have been called.
     *  @returns True only in case when the Step has caused the track to arrive
     *           at a boundary of a daughter. False, in all other cases.
     */
    inline G4bool EnteredDaughterVolume() const;

    /**
     * Verify if the step has exited the mother volume.
     */
    inline G4bool ExitedMotherVolume() const;

    /**
     * Run navigation in "check-mode", therefore using additional verifications
     * and more strict correctness conditions.
     *  @note Is effective only with G4VERBOSE set.
     */
    inline void CheckMode(G4bool mode);

    /**
     * Set/unset verbosity for pushed tracks (default is true).
     */
    inline G4bool IsCheckModeActive() const;
    inline void SetPushVerbosity(G4bool mode);

    /**
     * Prints the internal state of the Navigator (for debugging).
     * The level of detail is according to the verbosity.
     */
    void PrintState() const;

    /**
     * Obtains the transformations Global/Local (and inverse).
     *  @note Clients of these methods must copy the data
     *        if they need to keep it.
     */
    inline const G4AffineTransform& GetGlobalToLocalTransform() const;
    inline const G4AffineTransform  GetLocalToGlobalTransform() const;

    /**
     * Obtains mother to daughter transformation.
     */
    G4AffineTransform GetMotherToDaughterTransform(G4VPhysicalVolume* dVolume, 
                                                   G4int dReplicaNo,
                                                   EVolume dVolumeType );

    /**
     * Resets stack and minimum or navigator state machine necessary for reset
     * as needed by LocalGlobalPointAndSetup().
     *  @note Does not perform clears, resizes, or reset fLastLocatedPointLocal.
     */
    inline void ResetStackAndState();

    /**
     * Reports on severity of error and number of zero steps,
     * in case Navigator is stuck and is returning zero steps.
     * Values: 1 (small problem),  5 (correcting), 
     *         9 (ready to abandon), 10 (abandoned)
     *  @param[in,out] noZeroSteps Returns the number of zero steps in case
     *                 pointer is not null.
     *  @returns The error severity.
     */
    inline G4int SeverityOfZeroStepping( G4int* noZeroSteps ) const; 

    /**
     * Returns the local coordinate of the point in the reference system
     * of its containing volume that was found by LocalGlobalPointAndSetup().
     * The local coordinate of the last located track.
     */
    inline G4ThreeVector GetCurrentLocalCoordinate() const;

    /**
     * Computes and returns the local->global translation/rotation
     * of current volume.
     */
    inline G4ThreeVector NetTranslation() const;
    inline G4RotationMatrix NetRotation() const;

    /**
     * Enables best-possible evaluation of isotropic safety.
     */
    inline void EnableBestSafety( G4bool value = false );

    /**
     * Accessor & modifier for custom external navigation.
     */
    inline G4VExternalNavigation* GetExternalNavigation() const;
    void SetExternalNavigation(G4VExternalNavigation* externalNav);
   
    /**
     * Gets/sets alternative navigator for voxel volumes.
     */
    inline G4VoxelNavigation& GetVoxelNavigator();
    void SetVoxelNavigation(G4VoxelNavigation* voxelNav);

    /**
     * Cloning feature for use in MT applications to clone the navigator,
     * including external sub-navigator.
     *  @note Client has responsibility for ownership of the returned
     *        allocated pointer.
     *  @returns A pointer to the cloned navigator object.
     */
    inline G4Navigator* Clone() const;

    /**
     * Gets endpoint of last step.
     */
    inline G4ThreeVector GetLastStepEndPoint() const;

    /**
     * Derived navigators which rely on LocateGlobalPointAndSetup() need to
     * inform size of step, to maintain logic about arriving on boundary
     * for challenging cases.
     * Required in order to cope with multiple trials at boundaries
     * => Locate with use direction rather than simple, fast logic.
     */
    void InformLastStep(G4double lastStep,
                        G4bool entersDaughtVol,
                        G4bool exitsMotherVol );

  protected:

    /**
     * Saves the state: fValidExitNormal, fExitNormal, fExiting, fEntering,
     * fBlockedPhysicalVolume, fBlockedReplicaNo, fLastStepWasZero,
     * fLastLocatedPointLocal, fLocatedOutsideWorld, fEnteredDaughter,
     * fExitedMother, fPreviousSftOrigin, fPreviousSafety.
     */
    void SetSavedState();

    /**
     * Copy aspects of the state, to enable a non-state changing
     * call to ComputeStep().
     */
    void RestoreSavedState();
  
    /**
     * Utility method to reset the navigator state machine.
     */
    virtual void ResetState();

    /**
     * Returns position vector in local coordinate system, given a position
     * vector in world coordinate system.
     */
    inline G4ThreeVector ComputeLocalPoint(const G4ThreeVector& rGlobP) const;

    /**
     * Returns the local direction of the specified vector in the reference
     * system of the volume that was found by LocalGlobalPointAndSetup().
     * The Local Coordinates of point in world coordinate system.
     */
    inline G4ThreeVector ComputeLocalAxis(const G4ThreeVector& pVec) const;

    /**
     * Characterises the type of volume - normal/replicated/parameterised.
     */
    inline EVolume VolumeType(const G4VPhysicalVolume *pVol) const;

    /**
     * Characterises the daughters of given logical volume.
     */
    inline EVolume CharacteriseDaughters(const G4LogicalVolume *pLog) const;

    /**
     * Gets regular structure ID of first daughter.
     */
    inline G4int GetDaughtersRegularStructureId(const G4LogicalVolume *pLv) const;

    /**
     * Renavigates & resets hierarchy described by the current history:
     * Resets volumes and recomputes transforms and/or solids of
     * replicated/parameterised volumes.
     */
    virtual void SetupHierarchy();

    /**
     * Utility method to trigger overlaps check on a volume with reported
     * overlaps ordered by relevance. Used in ComputeStep() when loopings
     * with zero step are detected.
     */
    G4bool CheckOverlapsIterative(G4VPhysicalVolume* vol);

  private:

    /**
     * Logs and checks for steps larger than the tolerance.
     */
    void ComputeStepLog(const G4ThreeVector& pGlobalpoint,
                              G4double moveLenSq) const;

  protected:

    G4double kCarTolerance, fMinStep, fSqTol;  // Cached tolerances.

    // BEGIN State information ------------------------------------------------
    //

    /** Transformation and history of the current path
        through the geometrical hierarchy. */
    G4NavigationHistory fHistory;

    /** Endpoint of last ComputeStep().
        Can be used for optimisation (e.g. when computing safety). */
    G4ThreeVector fStepEndPoint;

    /** Position of the end-point of the last call to ComputeStep()
        in last local coordinates. */
    G4ThreeVector fLastStepEndPointLocal; 

    /** Verbosity level  [if > 0, printout can occur]. */
    G4int fVerbose = 0;
   
    /** A memory of whether in this Step a daughter volume is entered 
        (set in Compute & Locate).
        After Compute: it expects to enter a daughter
        After Locate:  it has entered a daughter. */
    G4bool fEnteredDaughter;

    /** A similar memory whether the Step exited current "mother" volume
        completely, not entering daughter. */
    G4bool fExitedMother;

    /** Set true if last Step was limited by geometry. */
    G4bool fWasLimitedByGeometry = false;

  private:

    /** Position of the last located point relative to its containing volume.
        This is coupled with the Boolean member 'fLocatedOutsideWorld'. */
    G4ThreeVector fLastLocatedPointLocal;

    /** Leaving volume normal, in the volume containing the exited
        volume's coordinate system.
        This is closely coupled with 'fValidExitNormal', which signals whether
        we have a (valid) normal for volume we're leaving. */
    G4ThreeVector fExitNormal;
   
    /** Leaving volume normal, in its own coordinate system. */
    G4ThreeVector fGrandMotherExitNormal;

    /** Leaving volume normal, in the global coordinate system. */
    G4ThreeVector fExitNormalGlobalFrame;

    /** Memory of last safety origin & value. Used in ComputeStep() to ensure
        that origin of current Step is in the same volume as the point of the
        last relocation. */
    G4ThreeVector fPreviousSftOrigin;
    G4double fPreviousSafety; 

    /** Memory of the mother volume during previous step.
        Intended use: inform user in case of stuck track. */
    G4VPhysicalVolume* fLastMotherPhys = nullptr;
   
    /** Identifies the volume and copy / replica number that is blocked
        (after exiting -- because the exit direction is along the exit normal)
        or a candidate for entry (after compute step). */
    G4VPhysicalVolume* fBlockedPhysicalVolume;
    G4int fBlockedReplicaNo;
   
    /** Count zero steps, as one or two can occur due to changing momentum at
        a boundary or at an edge common between volumes; several zero steps
        are likely a problem in the geometry description or in the navigation.
        Number of preceding moves that were Zero.
        Reset to 0 after finite step. */
    G4int fNumberZeroSteps;

    /** After this many failed/zero steps, act (push etc). */
    G4int fActionThreshold_NoZeroSteps = 10;  

    /** After this many failed/zero steps, abandon track. */
    G4int fAbandonThreshold_NoZeroSteps = 25; 

    /** States if the navigator is activated or not. */
    G4bool fActive = false;

    /** Whether ComputeStep() was called since the last call to a Locate().
        Uses: distinguish parts of state which differ before/after calls
        to ComputeStep() or one of the Locate() methods; avoid two consecutive
        calls to compute-step (illegal). */
    G4bool fLastTriedStepComputation = false; 

    /** Entering/Exiting volumes blocking/setup.
        o If exiting, volume ptr & replica number (set & used by Locate..())
          used for blocking on redescent of geometry;
        o If entering, volume ptr & replica number (set by ComputeStep(),
          used by Locate..()) of volume for 'automatic' entry. */
    G4bool fEntering, fExiting;
   
    /** Set true if have leaving volume normal. */
    G4bool fValidExitNormal;

    /** Whether the last ComputeStep moved Zero. Used to check for edges. */
    G4bool fLastStepWasZero;

    /** Whether the Navigator has detected an edge. */
    G4bool fLocatedOnEdge;       

    /** Whether the last call to Locate methods left the world. */
    G4bool fLocatedOutsideWorld;
   
    /** Whether frame is changed. */
    G4bool  fChangedGrandMotherRefFrame;

    /** Has it been computed since the last call to ComputeStep().
        Covers both Global and GrandMother. */
    G4bool  fCalculatedExitNormal;

    //
    // END State information --------------------------------------------------

    // Optional State information (created/used as needed)
    // 
   
    // Save key state information (NOT the navigation history stack)
    //
    struct G4SaveNavigatorState
    { 
       G4ThreeVector sExitNormal;  
       G4bool sValidExitNormal;    
       G4bool sEntering, sExiting;
       G4VPhysicalVolume* spBlockedPhysicalVolume;
       G4int sBlockedReplicaNo;  
       G4int sLastStepWasZero; 
       G4bool sWasLimitedByGeometry;

       //  Potentially relevant
       //
       G4bool sLocatedOutsideWorld;
       G4ThreeVector sLastLocatedPointLocal; 
       G4bool sEnteredDaughter, sExitedMother;
       G4ThreeVector sPreviousSftOrigin;
       G4double sPreviousSafety; 
    } fSaveState; 

    // ========================================================================
    // BEGIN -- Tracking Invariants

    /** A link to the topmost physical volume in the detector.
        Must be positioned at the origin and unrotated. */
    G4VPhysicalVolume* fTopPhysical = nullptr;

    // Helpers/Utility classes

    G4NormalNavigation fnormalNav;
    G4VoxelNavigation* fpvoxelNav;
    G4ParameterisedNavigation fparamNav;
    G4ReplicaNavigation freplicaNav;
    G4RegularNavigation fregularNav;
    G4VExternalNavigation* fpExternalNav = nullptr;
    G4VoxelSafety* fpVoxelSafety;
    G4SafetyCalculator* fpSafetyCalculator = nullptr;

    // Utility information

    /** Check-mode flag  [if true, more strict checks are performed]. */
    G4bool fCheck = false;

    /** Push flags  [if true, means a stuck particle has been pushed]. */
    G4bool fPushed = false, fWarnPush = true;

    // End -- Tracking Invariants
    // ========================================================================
};

#include "G4Navigator.icc"

#endif


// NOTES:
//
// The following methods provide detailed information when a Step has
// arrived at a geometrical boundary.  They distinguish between the different
// causes that can result in the track leaving its current volume.
//
// Four cases are possible:
//
// 1) The particle has reached a boundary of a daughter of the current volume:
//     (this could cause the relocation to enter the daughter itself
//     or a potential granddaughter or further descendant)
//     
// 2) The particle has reached a boundary of the current
//     volume, exiting into a mother (regardless the level
//     at which it is located in the tree):
//
// 3) The particle has reached a boundary of the current
//     volume, exiting into a volume which is not in its
//     parental hierarchy:
//
// 4) The particle is not on a boundary between volumes:
//     the function returns an exception, and the caller is
//     reccomended to compare the G4touchables associated
//     to the preStepPoint and postStepPoint to handle this case.
//
//   G4bool        EnteredDaughterVolume()
//   G4bool        IsExitNormalValid()
//   G4ThreeVector GetLocalExitNormal()
//
// The expected usefulness of these methods is to allow the caller to
// determine how to compute the surface normal at the volume boundary. The two
// possibilities are to obtain the normal from:
//
//   i) the solid associated with the volume of the initial point of the Step.
//      This is valid for cases 2 and 3.  
//      (Note that the initial point is generally the PreStepPoint of a Step).
//   or
// 
//  ii) the solid of the final point, ie of the volume after the relocation.
//      This is valid for case 1.
//      (Note that the final point is generally the PreStepPoint of a Step).
//
// This way the caller can always get a valid normal, pointing outside
// the solid for which it is computed, that can be used at his own
// discretion.
