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
// $Id: G4ITNavigator1.hh 85244 2014-10-27 08:24:13Z gcosmo $
//
// Original author: Paul Kent, July 95/96
//
/// \brief { Class description:
///
/// G4ITNavigator1 is a duplicate version of G4Navigator started from Geant4.9.5
/// initially written by Paul Kent and colleagues.
/// The only difference resides in the way the information is saved and managed
///
/// A class for use by the tracking management, able to obtain/calculate
/// dynamic tracking time information such as the distance to the next volume,
/// or to find the physical volume containing a given point in the world
/// reference system. The navigator maintains a transformation history and
/// other information to optimise the tracking time performance.}
//
// Contact : Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr)
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// We would be very happy hearing from you, send us your feedback! :)
//
// History:
// - Created.                                  Paul Kent,     Jul 95/96
// - Zero step protections                     J.A. / G.C.,   Nov  2004
// - Added check mode                          G. Cosmo,      Mar  2004
// - Made Navigator Abstract                   G. Cosmo,      Nov  2003
// - G4ITNavigator1 created                     M.K.,          Nov  2012
// *********************************************************************

#ifndef G4ITNAVIGATOR_HH
#define G4ITNAVIGATOR_HH

#include "geomdefs.hh"

#include "G4ThreeVector.hh"
#include "G4AffineTransform.hh"
#include "G4RotationMatrix.hh"

#include "G4LogicalVolume.hh"             // Used in inline methods
#include "G4GRSVolume.hh"                 //    "         "
#include "G4GRSSolid.hh"                  //    "         "
#include "G4TouchableHandle.hh"           //    "         "
#include "G4TouchableHistoryHandle.hh"

#include "G4NavigationHistory.hh"
#include "G4NormalNavigation.hh"
#include "G4VoxelNavigation.hh"
#include "G4ParameterisedNavigation.hh"
#include "G4ReplicaNavigation.hh"
#include "G4RegularNavigation.hh"

#include <iostream>

class G4VPhysicalVolume;


struct G4ITNavigatorState_Lock1
{
    virtual ~G4ITNavigatorState_Lock1(){;}
protected:
    G4ITNavigatorState_Lock1(){;}
};

class G4ITNavigator1
{
public:
	static const G4int fMaxNav = 8;   // rename to kMaxNoNav ??

  public:  // with description

  friend std::ostream& operator << (std::ostream &os, const G4ITNavigator1 &n);

  G4ITNavigator1();
    // Constructor - initialisers and setup.

  virtual ~G4ITNavigator1();
    // Destructor. No actions.

  // !>
    G4ITNavigatorState_Lock1* GetNavigatorState();
    void SetNavigatorState(G4ITNavigatorState_Lock1*);
    void NewNavigatorState();
  // <!

  virtual G4double ComputeStep(const G4ThreeVector &pGlobalPoint,
                               const G4ThreeVector &pDirection,
                               const G4double pCurrentProposedStepLength,
                                     G4double  &pNewSafety);
    // Calculate the distance to the next boundary intersected
    // along the specified NORMALISED vector direction and
    // from the specified point in the global coordinate
    // system. LocateGlobalPointAndSetup or LocateGlobalPointWithinVolume 
    // must have been called with the same global point prior to this call.
    // The isotropic distance to the nearest boundary is also
    // calculated (usually an underestimate). The current
    // proposed Step length is used to avoid intersection
    // calculations: if it can be determined that the nearest
    // boundary is >pCurrentProposedStepLength away, kInfinity
    // is returned together with the computed isotropic safety
    // distance. Geometry must be closed.

  G4double CheckNextStep(const G4ThreeVector &pGlobalPoint,
                         const G4ThreeVector &pDirection,
                         const G4double pCurrentProposedStepLength,
                               G4double &pNewSafety); 
    // Same as above, but do not disturb the state of the Navigator.

  virtual
  G4VPhysicalVolume* ResetHierarchyAndLocate(const G4ThreeVector &point,
                                             const G4ThreeVector &direction,
                                             const G4TouchableHistory &h);

    // Resets the geometrical hierarchy and search for the volumes deepest
    // in the hierarchy containing the point in the global coordinate space.
    // The direction is used to check if a volume is entered.
    // The search begin is the geometrical hierarchy at the location of the
    // last located point, or the endpoint of the previous Step if
    // SetGeometricallyLimitedStep() has been called immediately before.
    // 
    // Important Note: In order to call this the geometry MUST be closed.

  virtual
  G4VPhysicalVolume* LocateGlobalPointAndSetup(const G4ThreeVector& point,
                                             const G4ThreeVector* direction=0,
                                             const G4bool pRelativeSearch=true,
                                             const G4bool ignoreDirection=true);
    // Search the geometrical hierarchy for the volumes deepest in the hierarchy
    // containing the point in the global coordinate space. Two main cases are:
    //  i) If pRelativeSearch=false it makes use of no previous/state
    //     information. Returns the physical volume containing the point, 
    //     with all previous mothers correctly set up.
    // ii) If pRelativeSearch is set to true, the search begin is the
    //     geometrical hierarchy at the location of the last located point,
    //     or the endpoint of the previous Step if SetGeometricallyLimitedStep()
    //     has been called immediately before.
    // The direction is used (to check if a volume is entered) if either
    //   - the argument ignoreDirection is false, or
    //   - the Navigator has determined that it is on an edge shared by two or
    //     more volumes.  (This is state information.)
    // 
    // Important Note: In order to call this the geometry MUST be closed.

  virtual
  void LocateGlobalPointWithinVolume(const G4ThreeVector& position);
    // Notify the Navigator that a track has moved to the new Global point
    // 'position', that is known to be within the current safety.
    // No check is performed to ensure that it is within  the volume. 
    // This method can be called instead of LocateGlobalPointAndSetup ONLY if
    // the caller is certain that the new global point (position) is inside the
    // same volume as the previous position.  Usually this can be guaranteed
    // only if the point is within safety.

  inline void LocateGlobalPointAndUpdateTouchableHandle(
                const G4ThreeVector&       position,
                const G4ThreeVector&       direction,
                      G4TouchableHandle&   oldTouchableToUpdate,
                const G4bool               RelativeSearch = true);
    // First, search the geometrical hierarchy like the above method
    // LocateGlobalPointAndSetup(). Then use the volume found and its
    // navigation history to update the touchable.

  inline void LocateGlobalPointAndUpdateTouchable(
                const G4ThreeVector&       position,
                const G4ThreeVector&       direction,
                      G4VTouchable*        touchableToUpdate,
                const G4bool               RelativeSearch = true);
    // First, search the geometrical hierarchy like the above method
    // LocateGlobalPointAndSetup(). Then use the volume found and its
    // navigation history to update the touchable.

  inline void LocateGlobalPointAndUpdateTouchable(
                const G4ThreeVector&       position,
                      G4VTouchable*        touchableToUpdate,
                const G4bool               RelativeSearch = true);
    // Same as the method above but missing direction.

  inline void SetGeometricallyLimitedStep();
    // Inform the navigator that the previous Step calculated
    // by the geometry was taken in its entirety.

  virtual G4double ComputeSafety(const G4ThreeVector &globalpoint,
                                 const G4double pProposedMaxLength = DBL_MAX,
                                 const G4bool keepState = true);
    // Calculate the isotropic distance to the nearest boundary from the
    // specified point in the global coordinate system. 
    // The globalpoint utilised must be within the current volume.
    // The value returned is usually an underestimate.  
    // The proposed maximum length is used to avoid volume safety
    //  calculations.  The geometry must be closed.
    // To ensure minimum side effects from the call, keepState
    //  must be true.
  
  inline G4VPhysicalVolume* GetWorldVolume() const;
    // Return the current  world (`topmost') volume.

  inline void SetWorldVolume(G4VPhysicalVolume* pWorld);
    // Set the world (`topmost') volume. This must be positioned at
    // origin (0,0,0) and unrotated.

  inline G4GRSVolume* CreateGRSVolume() const;
  inline G4GRSSolid* CreateGRSSolid() const; 
  inline G4TouchableHistory* CreateTouchableHistory() const;
  inline G4TouchableHistory* CreateTouchableHistory(const G4NavigationHistory*) const;
    // `Touchable' creation methods: caller has deletion responsibility.

  virtual G4TouchableHistoryHandle CreateTouchableHistoryHandle() const;
    // Returns a reference counted handle to a touchable history.

  virtual G4ThreeVector GetLocalExitNormal(G4bool* valid);
  virtual G4ThreeVector GetLocalExitNormalAndCheck(const G4ThreeVector& point,
                                                         G4bool* valid);
  virtual G4ThreeVector GetGlobalExitNormal(const G4ThreeVector& point,
                                                  G4bool* valid);
    // Return Exit Surface Normal and validity too.
    // Can only be called if the Navigator's last Step has crossed a
    // volume geometrical boundary.
    // It returns the Normal to the surface pointing out of the volume that
    // was left behind and/or into the volume that was entered.
    // Convention:
    //   The *local* normal is in the coordinate system of the *final* volume.
    // Restriction:
    //   Normals are not available for replica volumes (returns valid= false)
    // These methods takes full care about how to calculate this normal,
    // but if the surfaces are not convex it will return valid=false.

  inline G4int GetVerboseLevel() const;
  inline void  SetVerboseLevel(G4int level);
    // Get/Set Verbose(ness) level.
    // [if level>0 && G4VERBOSE, printout can occur]

  inline G4bool IsActive() const;
    // Verify if the navigator is active.
  inline void  Activate(G4bool flag);
    // Activate/inactivate the navigator.

  inline G4bool EnteredDaughterVolume() const;
    // The purpose of this function is to inform the caller if the track is
    // entering a daughter volume while exiting from the current volume.
    // This method returns 
    // - True only in case 1) above, that is when the Step has caused
    //   the track to arrive at a boundary of a daughter.
    // - False in cases 2), 3) and 4), i.e. in all other cases.
    // This function is not guaranteed to work if SetGeometricallyLimitedStep()
    // was not called when it should have been called.
  inline G4bool ExitedMotherVolume() const;
    // Verify if the step has exited the mother volume.

  inline void   CheckMode(G4bool mode);
    // Run navigation in "check-mode", therefore using additional
    // verifications and more strict correctness conditions.
    // Is effective only with G4VERBOSE set.
  inline G4bool IsCheckModeActive() const;
  inline void   SetPushVerbosity(G4bool mode);
    // Set/unset verbosity for pushed tracks (default is true).

  void PrintState() const;
    // Print the internal state of the Navigator (for debugging).
    // The level of detail is according to the verbosity.

  inline const G4AffineTransform& GetGlobalToLocalTransform() const;
  inline const G4AffineTransform  GetLocalToGlobalTransform() const;
    // Obtain the transformations Global/Local (and inverse).
    // Clients of these methods must copy the data if they need to keep it.

  G4AffineTransform GetMotherToDaughterTransform(G4VPhysicalVolume* dVolume, 
                                                 G4int dReplicaNo,
                                                 EVolume dVolumeType );
    // Obtain mother to daughter transformation

  inline void ResetStackAndState();
    // Reset stack and minimum or navigator state machine necessary for reset
    // as needed by LocalGlobalPointAndSetup.
    // [Does not perform clears, resizes, or reset fLastLocatedPointLocal]

  inline G4int SeverityOfZeroStepping( G4int* noZeroSteps ) const; 
    // Report on severity of error and number of zero steps,
    // in case Navigator is stuck and is returning zero steps.
    // Values: 1 (small problem),  5 (correcting), 
    //         9 (ready to abandon), 10 (abandoned)

  void SetSavedState(); 
    // ( fValidExitNormal, fExitNormal, fExiting, fEntering, 
    //   fBlockedPhysicalVolume, fBlockedReplicaNo, fLastStepWasZero); 
  void RestoreSavedState(); 
    // Copy aspects of the state, to enable a non-state changing
    //  call to ComputeStep

  inline G4ThreeVector GetCurrentLocalCoordinate() const;
    // Return the local coordinate of the point in the reference system
    // of its containing volume that was found by LocalGlobalPointAndSetup.
    // The local coordinate of the last located track.

  inline G4ThreeVector NetTranslation() const;
  inline G4RotationMatrix NetRotation() const;
    // Compute+return the local->global translation/rotation of current volume.

  inline void EnableBestSafety( G4bool value= false );
    // Enable best-possible evaluation of isotropic safety

  virtual void ResetState();
    // Utility method to reset the navigator state machine.

 protected:  // with description

  inline G4ThreeVector ComputeLocalPoint(const G4ThreeVector& rGlobPoint) const;
    // Return position vector in local coordinate system, given a position
    // vector in world coordinate system.

  inline G4ThreeVector ComputeLocalAxis(const G4ThreeVector& pVec) const;
    // Return the local direction of the specified vector in the reference
    // system of the volume that was found by LocalGlobalPointAndSetup.
    // The Local Coordinates of point in world coordinate system.

  inline EVolume VolumeType(const G4VPhysicalVolume *pVol) const;
    // Characterise `type' of volume - normal/replicated/parameterised.

  inline EVolume CharacteriseDaughters(const G4LogicalVolume *pLog) const;
    // Characterise daughter of logical volume.

  inline G4int GetDaughtersRegularStructureId(const G4LogicalVolume *pLog) const;
    // Get regular structure ID of first daughter

  virtual void SetupHierarchy();
    // Renavigate & reset hierarchy described by current history
    // o Reset volumes
    // o Recompute transforms and/or solids of replicated/parameterised
    //   volumes.

 private:

  G4ITNavigator1(const G4ITNavigator1&);
  G4ITNavigator1& operator=(const G4ITNavigator1&);
    // Private copy-constructor and assignment operator.

  void ComputeStepLog(const G4ThreeVector& pGlobalpoint,
                            G4double moveLenSq) const;
    // Log and checks for steps larger than the tolerance

 protected:  // without description

  G4double kCarTolerance;
    // Geometrical tolerance for surface thickness of shapes.

  //
  // BEGIN State information
  //

  G4NavigationHistory fHistory;
    // Transformation and history of the current path
    // through the geometrical hierarchy.

  G4bool fEnteredDaughter;
    // A memory of whether in this Step a daughter volume is entered 
    // (set in Compute & Locate).
    //  After Compute: it expects to enter a daughter
    //  After Locate:  it has entered a daughter

  G4bool fExitedMother;
    // A similar memory whether the Step exited current "mother" volume
    // completely, not entering daughter.

  G4bool fWasLimitedByGeometry;
    // Set true if last Step was limited by geometry.

  G4ThreeVector fStepEndPoint;
    // Endpoint of last ComputeStep 
    // can be used for optimisation (e.g. when computing safety).
  G4ThreeVector fLastStepEndPointLocal; 
    // Position of the end-point of the last call to ComputeStep 
    // in last Local coordinates.

  G4int  fVerbose;
    // Verbose(ness) level  [if > 0, printout can occur].

 private:

  G4bool fActive;
    // States if the navigator is activated or not.

  G4bool fLastTriedStepComputation; 
    // Whether ComputeStep was called since the last call to a Locate method
    // Uses: - distinguish parts of state which differ before/after calls
    //         to ComputeStep or one of the Locate methods;
    //       - avoid two consecutive calls to compute-step (illegal).

  G4bool fEntering,fExiting;
    // Entering/Exiting volumes blocking/setup
    // o If exiting
    //      volume ptr & replica number (set & used by Locate..())
    //      used for blocking on redescent of geometry
    // o If entering
    //      volume ptr & replica number (set by ComputeStep(),used by
    //      Locate..()) of volume for `automatic' entry

  G4VPhysicalVolume *fBlockedPhysicalVolume;
  G4int fBlockedReplicaNo;

  G4ThreeVector fLastLocatedPointLocal;
    // Position of the last located point relative to its containing volume.
  G4bool fLocatedOutsideWorld;
    // Whether the last call to Locate methods left the world

  G4bool fValidExitNormal;    // Set true if have leaving volume normal
  G4ThreeVector fExitNormal;  // Leaving volume normal, in the
                              // volume containing the exited
                              // volume's coordinate system
  G4ThreeVector fGrandMotherExitNormal;  // Leaving volume normal, in its 
                                         // own coordinate system
  G4bool  fChangedGrandMotherRefFrame;   // Whether frame is changed

  G4ThreeVector fExitNormalGlobalFrame;  // Leaving volume normal, in the
                                         // global coordinate system
  G4bool  fCalculatedExitNormal;  // Has it been computed since
                                  // the last call to ComputeStep
                                  // Covers both Global and GrandMother
   
  // Count zero steps - as one or two can occur due to changing momentum at
  //                    a boundary or at an edge common between volumes
  //                  - several are likely a problem in the geometry
  //                    description or in the navigation
  //
  G4bool fLastStepWasZero;
    // Whether the last ComputeStep moved Zero. Used to check for edges.

  G4bool fLocatedOnEdge;       
    // Whether the Navigator has detected an edge
  G4int fNumberZeroSteps;
    // Number of preceding moves that were Zero. Reset to 0 after finite step
  G4int fActionThreshold_NoZeroSteps;  
    // After this many failed/zero steps, act (push etc) 
  G4int fAbandonThreshold_NoZeroSteps; 
    // After this many failed/zero steps, abandon track

  G4ThreeVector  fPreviousSftOrigin;
  G4double       fPreviousSafety; 
    // Memory of last safety origin & value. Used in ComputeStep to ensure
    // that origin of current Step is in the same volume as the point of the
    // last relocation

  //
  // END State information
  //

  // Save key state information (NOT the navigation history stack)
  //
  struct G4SaveNavigatorState : public G4ITNavigatorState_Lock1
  { 
     G4SaveNavigatorState();
     virtual ~G4SaveNavigatorState(){;}
     G4ThreeVector sExitNormal;  
     G4bool sValidExitNormal;    
     G4bool sEntering, sExiting;
     G4VPhysicalVolume* spBlockedPhysicalVolume;
     G4int sBlockedReplicaNo;  
     G4int sLastStepWasZero; 

     // !>
     G4bool sLocatedOnEdge;
     G4bool sWasLimitedByGeometry;
     G4bool sPushed;
     G4int  sNumberZeroSteps;
     // <!

     //  Potentially relevant
     //
     G4bool sLocatedOutsideWorld;
     G4ThreeVector sLastLocatedPointLocal; 
     G4bool sEnteredDaughter, sExitedMother;
     G4ThreeVector  sPreviousSftOrigin;
     G4double       sPreviousSafety; 
  } ;

  G4SaveNavigatorState* fpSaveState;


  // Tracking Invariants
  //
  G4VPhysicalVolume  *fTopPhysical;
    // A link to the topmost physical volume in the detector.
    // Must be positioned at the origin and unrotated.

  // Utility information
  //
  G4bool fCheck;
    // Check-mode flag  [if true, more strict checks are performed].
  G4bool fPushed, fWarnPush;
    // Push flags  [if true, means a stuck particle has been pushed].

  // Helpers/Utility classes
  //
  G4NormalNavigation  fnormalNav;
  G4VoxelNavigation fvoxelNav;
  G4ParameterisedNavigation fparamNav;
  G4ReplicaNavigation freplicaNav;
  G4RegularNavigation fregularNav;
  G4VoxelSafety       *fpVoxelSafety;
};

#include "G4ITNavigator1.icc"

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
