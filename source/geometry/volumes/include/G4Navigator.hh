// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Navigator.hh,v 1.4 1999-12-15 16:40:11 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4Navigator Paul Kent July 95/96
//
// A class for use by the tracking management, able to obtain/calculate
// dynamic tracking time information such as the distance to the next volume,
// or to find the physical volume containing a given point in the world
// reference system. The navigator maintains a transformation history and
// other information to optimise the tracking time performance.
//
// Member functions:
//
// G4Navigator()
//    Constructor. No actions yet.
//
// ~G4Navigator()
//    Destructor. No actions yet.
//
// G4double ComputeStep(const G4ThreeVector &globalpoint,
//                      const G4ThreeVector &direction,
// 		        const G4double pCurrentProposedStepLength,
//		        G4double &newSafety)
// Calculate the distance to the next boundary intersected
// along the specified NORMALISED vector direction and
// from the specified point in the global coordiante
// system.  LocateGlobalPointAndSetup or LocateGlobalPointWithinVolume 
// must have been called with the same global point prior to this call.
// The isotropic distance to the nearest boundary is also
// calculated (usually an underestimate). The current
// proposed Step length is used to avoid intersection
// calculations: if it can be determined that the nearest
// boundary is >pCurrentProposedStepLength away, kInfinity
// is returned together with the computed isotropic safety
// distance. Geometry must be closed.
//
// G4double ComputeSafety(const G4ThreeVector &globalpoint,
// 		          const G4double pProposedMaxLength=DBL_MAX )
//
//  Calculate the isotropic distance to the nearest boundary from the
// specified point in the global coordinate system. 
// The globalpoint utilised must be within the current volume.
// The value returned is usually an underestimate.  
// The proposed maximum length is used to avoid volume safety
// calculations.  The geometry must be closed.
//
// G4ThreeVector GetCurrentLocalCoordinate() const
//   Return the local coordinate of the point in the reference system
//   of its containing volume that was found by
//   LocalGlobalPointAndSetup.
//
// G4ThreeVector ComputeLocalAxis(const G4ThreeVector& pGVector) const
//   Return the local direction of the specified vector in the reference
//   system of the volume that was found by LocalGlobalPointAndSetup.
//
// G4VPhysicalVolume* GetWorldVolume() const
//   Return the current  world (`topmost') volume
//
// G4VPhysicalVolume* LocateGlobalPointAndSetup(const G4ThreeVector& point,
//            				        const G4ThreeVector* direction=0,
// 					       G4bool pRelativeSearch=true)

// Search the geometrical hierarchy for the volumes deepest in the hierarchy
// containing the point in the global coordinate space.  Two main cases are:
//  i) If pRelativeSearch=false it makes use of no previous/state
//     information. Returns the physical volume containing the point, 
//     with all previous mothers correctly set up.
// ii) If pRelativeSearch is set to true, the search begin is the geometrical
//     hierarchy at the location of the last located point, or the endpoint of
//     the previous Step if SetGeometricallyLimitedStep() has been called
//     immediately before.
// The direction is used (to check if a volume is entered) only if 
// the Navigator has determined that it is on an edge shared by two or more
// volumes.  (This is state information.)
// The geometry must be closed.
//
//
// void LocateGlobalPointWithinVolume( const  G4ThreeVector& position)
//   Notify the Navigator that a track has moved to the new Global point
//   'position', that is know to be within the current safety.
//   No check is performed to ensure that it is within  the volume. 
//   This method can be called instead of LocateGlobalPointAndSetup  ONLY if
//   the caller is certain that the new global point (position) is inside the
//   same volume as the previous position.  Usually this can be guaranteed
//   only if the point is within safety.
//
// void   
// LocateGlobalPointAndUpdateTouchable( const G4ThreeVector& position,
//				        const G4ThreeVector& globalDirection,
//                                      G4Touchable*    touchableToUpdate,
//                                      const G4bool    relativeSearch=true);
// First, search the geometrical hierarchy like the above method
// LocateGlobalPointAndSetup(const G4ThreeVector&, G4bool).
// Then use the volume found and its navigation history to
// update the touchable.
//

// void SetGeometricallyLimitedStep()
//   Inform the navigator that the previous Step calculated
//   by the geometry was taken in its entirety
//
// void SetWorldVolume(G4VPhysicalVolume* pWorld)
//   Set the world (`topmost') volume. This must be
//   positioned at (0,0,0) and unrotated.
//
// G4ThreeVector GetLocalExitNormal(G4bool* valid);
//     Can be called only if the Navigator's last Step has arrived at 
//      a geometrical boundary.  
//     It returns the Normal to the surface pointing out of the volume that
//      was left behind and/or into the volume that was entered.
//      (The normal is in the coordinate system of the final volume.)
//     This function takes full care about how to calculate this normal,
//      but if the surfaces are not convex it will return valid=false.
//
// The following methods provide (more detailed) information when a Step has
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
//  G4bool        EnteredDaughterVolume();
//   The purpose of this function is to inform the caller if the track is
//   entering a daughter volume while exiting from the current volume.
//    This method returns 
//     - True only in case 1) above, that is when
//     the Step has caused the track to arrive at a boundary of a daughter.
//     - False in cases 2), 3) and 4), ie in all other cases.
//   This function is not guaranteed to work if SetGeometricallyLimitedStep()
//    was not called when it should have been called.
//
//  G4bool        IsExitNormalValid();
//    Return true if the Navigator
//     i) found that it is exiting the previous volume and is not entering a
//        daughter volume.
//    ii) has obtained the normal for the previous volume.
//
//  G4ThreeVector GetLocalExitNormal(); 
//     Can be called (ie can return a valid result) only if Navigator replied 
//      true to IsExitNormalValid()
//     It returns the ExitNormal of the previous volume
//       (The normal is in the coordinate system of the final volume.)
//
// The expected usefulness of these methods is to allow the caller to
// determine how to compute the surface normal at the volume boundary. The two
// possibilities are to obtain the normal from:
//
//   i) the solid associated with the volume of the initial point of the Step.
//      This is valid for cases 2 and 3.  
//        (Note that the initial point is generally the PreStepPoint of a Step).
//   or 
//  ii) the solid of the final point, ie of the volume after the relocation.
//      This is valid for case 1.
//        (Note that the final point is generally the PreStepPoint of a Step).
//
// This way the caller can always get a valid normal, pointing outside
// the solid for which it is computed, that can be used at his own
// discretion.
//

// 18th Dec 1997  John Allison  Temporary change -> NON-CONST (8 functions)

#ifndef G4NAVIGATOR_HH
#define G4NAVIGATOR_HH
#include "geomdefs.hh"

#include "G4ThreeVector.hh"
#include "G4AffineTransform.hh"
#include "G4RotationMatrix.hh"

#include "G4VPhysicalVolume.hh"
#include "G4VPVParameterisation.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"

#include "G4GRSVolume.hh"
#include "G4GRSSolid.hh"
#include "G4TouchableHistory.hh"

#include "G4NavigationHistory.hh"
#include "G4NormalNavigation.hh"
#include "G4VoxelNavigation.hh"
#include "G4ParameterisedNavigation.hh"
#include "G4ReplicaNavigation.hh"

#include "g4std/iostream"

class G4Navigator
{
public:
  friend G4std::ostream& operator << (G4std::ostream &os, const G4Navigator &n);

// Constructor - initialisers and setup
  G4Navigator();
  ~G4Navigator();

// Compute the next geometric Step
  G4double ComputeStep(const G4ThreeVector &pGlobalPoint,
		       const G4ThreeVector &pDirection,
		       const G4double pCurrentProposedStepLength,
		       G4double	&pNewSafety);
    
// Locate the point in the hierarchy
  G4VPhysicalVolume* LocateGlobalPointAndSetup(const G4ThreeVector& point,
					       const G4ThreeVector* direction=0,
					       const G4bool pRelativeSearch=true);

  G4VPhysicalVolume* LocateGlobalPointAndSetup(const G4ThreeVector &p,
  			   const G4TouchableHistory &h);

  void LocateGlobalPointWithinVolume( const  G4ThreeVector& position);

  void               LocateGlobalPointAndUpdateTouchable(
			  const G4ThreeVector&       position,
			  const G4ThreeVector&       direction,
				G4VTouchable*        touchableToUpdate,
			  const G4bool               RelativeSearch  =true);

  // Old version (missing direction) : not recommended
  //  replace with new version above.

  void               LocateGlobalPointAndUpdateTouchable(
			  const G4ThreeVector&       position,
				G4VTouchable*        touchableToUpdate,
			  const G4bool               RelativeSearch  =true);

// Inform the navigator that the previous Step calculated
// by the geometry was taken in its entirety
  void SetGeometricallyLimitedStep();

//  Calculate the isotropic distance to the nearest boundary from the
// specified point in the global coordinate system. 
G4double ComputeSafety(const G4ThreeVector &globalpoint,
 		       const G4double pProposedMaxLength=DBL_MAX );

// Return the local coordinate of the last located track
  G4ThreeVector GetCurrentLocalCoordinate() const;
// Return position vector in local coordinate system, given a position vector 
//  in world coord system
  G4ThreeVector ComputeLocalPoint(const G4ThreeVector& rGlobPoint) const;
// Return Local Coordinates of point in world coordinate system
//   
  G4ThreeVector ComputeLocalAxis(const G4ThreeVector& pVec) const;

// Compute+return the local->global translation/rotation of current volume
  G4ThreeVector NetTranslation() const;
  G4RotationMatrix NetRotation() const;

// Return the current  world (`topmost') volume
  G4VPhysicalVolume* GetWorldVolume() const;

// Set the world (`topmost') volume. Check at origin and unrotated.
  void SetWorldVolume(G4VPhysicalVolume* pWorld);

// `Touchable' creation methods: caller has deletion responsibility
  G4GRSVolume* CreateGRSVolume() const;
  G4GRSSolid* CreateGRSSolid() const; 
  G4TouchableHistory* CreateTouchableHistory() const;

//  Obtain the Normal vector to a surface (in local coordinates)
  G4bool        IsExitNormalValid();
  G4ThreeVector GetLocalExitNormal(); // Valid only if Navigator replied true
		        // It is in the coordinate system of the final volume

  // More powerful calculation of Exit Surface Normal, returns validity too.
  //  But it can only be called if the Step has crossed a volume boundary.
  G4ThreeVector GetLocalExitNormal(G4bool* valid);

  // 
  G4bool        EnteredDaughterVolume();
  // G4bool        ExitedMotherVolume();

  // Get/Set Verbose(ness) level (if level>0 && G4VERBOSE, printout can occur)
  G4int GetVerboseLevel();
  void  SetVerboseLevel(G4int level);

  //  Print the internal state of the Navigator (for debugging).
  //   The level of detail is according to the verbosity.
  void  PrintState();      

//  Obtain the transformations  Global/Local (and inverse)
//   Clients of these methods must copy the data if they need to keep it.
  const G4AffineTransform& GetGlobalToLocalTransform() const;
  const G4AffineTransform  GetLocalToGlobalTransform() const;

protected:
// Reset stack and minimum or navigator state machine necessary for reset
// as needed by LocalGlobalPointAnsSetup
// [Does not perform clears, resizes, or reset fLastLocatedPointLocal]
  void ResetStackAndState();

// Characterise `type' of volume - normal/replicated/parameterised
  EVolume VolumeType(const G4VPhysicalVolume *pVol) const;
// Characterise daughter of logical volume
  EVolume CharacteriseDaughters(const G4LogicalVolume *pLog) const;

// Renavigate & reset hierarchy described by current history
// o Reset volumes
// o Recompute transforms and/or solids of replicated/parameterised vols
  void SetupHierarchy();
private:

//
// BEGIN State information
//

// Position of the last located point relative to its containing volume
  G4ThreeVector fLastLocatedPointLocal;
  
// Set true if last Step was limited by geometry.
  G4bool fWasLimitedByGeometry;
  G4bool fEntering,fExiting;

// Entering/Exiting volumes blocking/setup

// o If exiting
//      volume ptr & replica number (set & used by Locate..())
//      used for blocking on redescent of geometry
  G4VPhysicalVolume *fBlockedPhysicalVolume;
  G4int fBlockedReplicaNo;

// o If entering
//      volume ptr & replica number (set by ComputeStep(),used by
//      Locate..()) of volume for `automatic' entry
  G4VPhysicalVolume *fCandidatePhysicalVolume;
  G4int fCandidateReplicaNo;

  G4bool fEnteredDaughter;    // A memory of whether in this Step a daughter
                              //  volume is entered (set in Compute & Locate)
			      //  After Compute: it expects to enter a daughter
                              //  After Locate:  it has entered a daughter
  G4bool fExitedMother;       // A similar memory whether the Step exited 
                              //  current "mother" volume completely, 
                              //  not entering daughter.
  
  G4bool fValidExitNormal;	// Set true if have leaving volume normal
  G4ThreeVector fExitNormal;	// Leaving volume normal, in the
				// volume containing the exited
				// volume's coordinate system

  G4NavigationHistory fHistory;
				// Transformation & `path' history
                                // of current path through geomtrical
                                // hierarchy

  G4bool  fLastStepWasZero;
                                // Whether the last ComputeStep moved Zero
                                //  Used to check for edges.
  G4bool  fLocatedOnEdge;       
                                // Whether the Navigator has detected an edge

  G4ThreeVector  fPreviousSftOrigin;
  G4double       fPreviousSafety; 
                                // Memory of last safety origin & value.
                                //  Used in ComputeStep to ensure that
                                //  origin of current Step is in the same
                                //  volume as the point of the last relocation.
//
// END State information
//

//
// BEGIN Tracking Invariants
//

// A link to the topmost physical volume in the detector.
// Must be positioned at the origin and unrotated.
  G4VPhysicalVolume	*fTopPhysical;
//
// END Tracking Invariants
//

//
// BEGIN Utility information
//
  G4int  fVerbose;   // Verbose(ness) level  (if > 0, printout can occur)
//
// END Utility Invariants
//

// 
// Helpers/Utility classes
//
  G4NormalNavigation  fnormalNav;
  G4VoxelNavigation fvoxelNav;
  G4ParameterisedNavigation fparamNav;
  G4ReplicaNavigation freplicaNav;
};

#include "G4Navigator.icc"

#endif






