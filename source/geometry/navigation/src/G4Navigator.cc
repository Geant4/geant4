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
//
// $Id: G4Navigator.cc 2012-11-01 japost Exp $
// GEANT4 tag $ Name:  $
// 
// class G4Navigator Implementation
//
// Original author: Paul Kent, July 95/96
// Responsible 1996-present: John Apostolakis
// Co-maintainer:            Gabriele Cosmo
// Additional revisions by: Pedro Arce, Vladimir Grichine
// --------------------------------------------------------------------

#include <iomanip>

#include "G4Navigator.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4GeometryTolerance.hh"
#include "G4VPhysicalVolume.hh"

#include "G4VoxelSafety.hh"

// ********************************************************************
// Constructor
// ********************************************************************
//
G4Navigator::G4Navigator()
  : fWasLimitedByGeometry(false), fVerbose(0),
    fTopPhysical(0), fCheck(false), fPushed(false), fWarnPush(true)
{
  fActive= false; 
  fLastTriedStepComputation= false;

  ResetStackAndState();
    // Initialises also all 
    // - exit / entry flags
    // - flags & variables for exit normals
    // - zero step counters
    // - blocked volume 

  fActionThreshold_NoZeroSteps  = 10; 
  fAbandonThreshold_NoZeroSteps = 25; 

  kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  fMinStep = 0.05*kCarTolerance;
  fSqTol = kCarTolerance*kCarTolerance;

  fregularNav.SetNormalNavigation( &fnormalNav );

  fStepEndPoint = G4ThreeVector( kInfinity, kInfinity, kInfinity ); 
  fLastStepEndPointLocal = G4ThreeVector( kInfinity, kInfinity, kInfinity ); 

  fpVoxelSafety= new G4VoxelSafety();
}

// ********************************************************************
// Destructor
// ********************************************************************
//
G4Navigator::~G4Navigator()
{
  delete fpVoxelSafety;
}

// ********************************************************************
// ResetHierarchyAndLocate
// ********************************************************************
//
G4VPhysicalVolume*
G4Navigator::ResetHierarchyAndLocate(const G4ThreeVector &p,
                                     const G4ThreeVector &direction,
                                     const G4TouchableHistory &h)
{
  ResetState();
  fHistory = *h.GetHistory();
  SetupHierarchy();
  fLastTriedStepComputation= false;  // Redundant, but best
  return LocateGlobalPointAndSetup(p, &direction, true, false);
}

// ********************************************************************
// LocateGlobalPointAndSetup
//
// Locate the point in the hierarchy return 0 if outside
// The direction is required 
//    - if on an edge shared by more than two surfaces 
//      (to resolve likely looping in tracking)
//    - at initial location of a particle
//      (to resolve potential ambiguity at boundary)
// 
// Flags on exit: (comments to be completed)
// fEntering         - True if entering `daughter' volume (or replica)
//                     whether daughter of last mother directly 
//                     or daughter of that volume's ancestor.
// fExiting          - True if exited 'mother' volume
//                     (always ? - how about if going back down ? - tbc)
// ********************************************************************
//
G4VPhysicalVolume* 
G4Navigator::LocateGlobalPointAndSetup( const G4ThreeVector& globalPoint,
                                        const G4ThreeVector* pGlobalDirection,
                                        const G4bool relativeSearch,
                                        const G4bool ignoreDirection )
{
  G4bool notKnownContained=true, noResult;
  G4VPhysicalVolume *targetPhysical;
  G4LogicalVolume *targetLogical;
  G4VSolid *targetSolid=0;
  G4ThreeVector localPoint, globalDirection;
  EInside insideCode;
  
  G4bool considerDirection = (!ignoreDirection) || fLocatedOnEdge;

  fLastTriedStepComputation=   false;   
  fChangedGrandMotherRefFrame= false;  // For local exit normal
   
  if( considerDirection && pGlobalDirection != 0 )
  {
    globalDirection=*pGlobalDirection;
  }

#ifdef G4VERBOSE
  if( fVerbose > 2 )
  {
    G4int oldcoutPrec = G4cout.precision(8);
    G4cout << "*** G4Navigator::LocateGlobalPointAndSetup: ***" << G4endl; 
    G4cout << "    Called with arguments: " << G4endl
           << "        Globalpoint = " << globalPoint << G4endl
           << "        RelativeSearch = " << relativeSearch  << G4endl;
    if( fVerbose >= 4 )
    {
      G4cout << "    ----- Upon entering:" << G4endl;
      PrintState();
    }
    G4cout.precision(oldcoutPrec);
  }
#endif

  G4int noLevelsExited=0 ;
  G4int noLevelsEntered= 0;

  if ( !relativeSearch )
  {
    ResetStackAndState();
  }
  else
  {
    if ( fWasLimitedByGeometry )
    {
      fWasLimitedByGeometry = false;
      fEnteredDaughter = fEntering;   // Remember
      fExitedMother = fExiting;       // Remember
      if ( fExiting )
      {
        noLevelsExited++;  // count this first level entered too

        if ( fHistory.GetDepth() )
        {
          fBlockedPhysicalVolume = fHistory.GetTopVolume();
          fBlockedReplicaNo = fHistory.GetTopReplicaNo();
          fHistory.BackLevel();
        }
        else
        {
          fLastLocatedPointLocal = localPoint;
          fLocatedOutsideWorld = true;
          fBlockedPhysicalVolume = 0;   // to be sure
          fBlockedReplicaNo = -1;
          fEntering = false;            // No longer 
          fEnteredDaughter = false; 
          fExitedMother = true;      // ??
          
          return 0;           // Have exited world volume
        }
        // A fix for the case where a volume is "entered" at an edge
        // and a coincident surface exists outside it.
        //  - This stops it from exiting further volumes and cycling
        //  - However ReplicaNavigator treats this case itself
        //
        // assert( fBlockedPhysicalVolume!=0 );

        // Expect to be on edge => on surface
        //
        if ( fLocatedOnEdge && (VolumeType(fBlockedPhysicalVolume)!=kReplica ))
        { 
          fExiting= false;
          // Consider effect on Exit Normal !?
        }
      }
      else
        if ( fEntering )
        {
          // assert( fBlockedPhysicalVolume!=0 );

          noLevelsEntered++;   // count the first level entered too

          switch (VolumeType(fBlockedPhysicalVolume))
          {
            case kNormal:
              fHistory.NewLevel(fBlockedPhysicalVolume, kNormal,
                                fBlockedPhysicalVolume->GetCopyNo());
              break;
            case kReplica:
              freplicaNav.ComputeTransformation(fBlockedReplicaNo,
                                                fBlockedPhysicalVolume);
              fHistory.NewLevel(fBlockedPhysicalVolume, kReplica,
                                fBlockedReplicaNo);
              fBlockedPhysicalVolume->SetCopyNo(fBlockedReplicaNo);
              break;
            case kParameterised:
              if( fBlockedPhysicalVolume->GetRegularStructureId() == 0 )
              {
                G4VSolid *pSolid;
                G4VPVParameterisation *pParam;
                G4TouchableHistory parentTouchable( fHistory );
                pParam = fBlockedPhysicalVolume->GetParameterisation();
                pSolid = pParam->ComputeSolid(fBlockedReplicaNo,
                                              fBlockedPhysicalVolume);
                pSolid->ComputeDimensions(pParam, fBlockedReplicaNo,
                                          fBlockedPhysicalVolume);
                pParam->ComputeTransformation(fBlockedReplicaNo,
                                              fBlockedPhysicalVolume);
                fHistory.NewLevel(fBlockedPhysicalVolume, kParameterised,
                                  fBlockedReplicaNo);
                fBlockedPhysicalVolume->SetCopyNo(fBlockedReplicaNo);
                //
                // Set the correct solid and material in Logical Volume
                //
                G4LogicalVolume *pLogical;
                pLogical = fBlockedPhysicalVolume->GetLogicalVolume();
                pLogical->SetSolid( pSolid );
                pLogical->UpdateMaterial(pParam ->
                  ComputeMaterial(fBlockedReplicaNo,
                                  fBlockedPhysicalVolume, 
                                  &parentTouchable));
              }
              break;
          }
          fEntering = false;
          fBlockedPhysicalVolume = 0;
          localPoint = fHistory.GetTopTransform().TransformPoint(globalPoint);
          notKnownContained = false;
        }
    }
    else
    {
      fBlockedPhysicalVolume = 0;
      fEntering = false;
      fEnteredDaughter = false;  // Full Step was not taken, did not enter
      fExiting = false;
      fExitedMother = false;     // Full Step was not taken, did not exit
    }
  }
  //
  // Search from top of history up through geometry until
  // containing volume found:
  // If on 
  // o OUTSIDE - Back up level, not/no longer exiting volumes
  // o SURFACE and EXITING - Back up level, setting new blocking no.s
  // else
  // o containing volume found
  //

  while (notKnownContained)  // Loop checking, 07.10.2016, J.Apostolakis
  {
    if ( fHistory.GetTopVolumeType()!=kReplica )
    {
      targetSolid = fHistory.GetTopVolume()->GetLogicalVolume()->GetSolid();
      localPoint = fHistory.GetTopTransform().TransformPoint(globalPoint);
      insideCode = targetSolid->Inside(localPoint);
#ifdef G4VERBOSE
      if(( fVerbose == 1 ) && ( fCheck ))
      {
         G4String solidResponse = "-kInside-";
         if (insideCode == kOutside)
           solidResponse = "-kOutside-";
         else if (insideCode == kSurface)
           solidResponse = "-kSurface-";
         G4cout << "*** G4Navigator::LocateGlobalPointAndSetup(): ***" << G4endl
                << "    Invoked Inside() for solid: " << targetSolid->GetName()
                << ". Solid replied: " << solidResponse << G4endl
                << "    For local point p: " << localPoint << G4endl;
      }
#endif
    }
    else
    {
      insideCode = freplicaNav.BackLocate(fHistory, globalPoint, localPoint,
                                          fExiting, notKnownContained);
      // !CARE! if notKnownContained returns false then the point is within
      // the containing placement volume of the replica(s). If insidecode
      // will result in the history being backed up one level, then the
      // local point returned is the point in the system of this new level
    }


    if ( insideCode==kOutside )
    {
      noLevelsExited++; 
      if ( fHistory.GetDepth() )
      {
        fBlockedPhysicalVolume = fHistory.GetTopVolume();
        fBlockedReplicaNo = fHistory.GetTopReplicaNo();
        fHistory.BackLevel();
        fExiting = false;

        if( noLevelsExited > 1 )
        {
          // The first transformation was done by the sub-navigator
          //
          const G4RotationMatrix* mRot = fBlockedPhysicalVolume->GetRotation();
          if( mRot )
          { 
            fGrandMotherExitNormal *= (*mRot).inverse();
            fChangedGrandMotherRefFrame= true;
          }
        }
      }
      else
      {
        fLastLocatedPointLocal = localPoint;
        fLocatedOutsideWorld = true;
          // No extra transformation for ExitNormal - is in frame of Top Volume
        return 0;         // Have exited world volume
      }
    }
    else
      if ( insideCode==kSurface )
      {
        G4bool isExiting = fExiting;
        if( (!fExiting)&&considerDirection )
        {
          // Figure out whether we are exiting this level's volume
          // by using the direction
          //
          G4bool directionExiting = false;
          G4ThreeVector localDirection =
              fHistory.GetTopTransform().TransformAxis(globalDirection);

          // Make sure localPoint in correct reference frame
          //     ( Was it already correct ? How ? )
          //
          localPoint= fHistory.GetTopTransform().TransformPoint(globalPoint);
          if ( fHistory.GetTopVolumeType()!=kReplica )
          {
            G4ThreeVector normal = targetSolid->SurfaceNormal(localPoint);
            directionExiting = normal.dot(localDirection) > 0.0;
            isExiting = isExiting || directionExiting;
          }
        }
        if( isExiting )
        {
          noLevelsExited++; 
          if ( fHistory.GetDepth() )
          {
            fBlockedPhysicalVolume = fHistory.GetTopVolume();
            fBlockedReplicaNo = fHistory.GetTopReplicaNo();
            fHistory.BackLevel();
            //
            // Still on surface but exited volume not necessarily convex
            //
            fValidExitNormal = false;

            if( noLevelsExited > 1 )
            {
              // The first transformation was done by the sub-navigator
              //
              const G4RotationMatrix* mRot =
                    fBlockedPhysicalVolume->GetRotation();
              if( mRot )
              { 
                fGrandMotherExitNormal *= (*mRot).inverse();
                fChangedGrandMotherRefFrame= true;
              }
            }
          } 
          else
          {
            fLastLocatedPointLocal = localPoint;
            fLocatedOutsideWorld = true;
              // No extra transformation for ExitNormal, is in frame of Top Vol
            return 0;          // Have exited world volume
          }
        }
        else
        {
          notKnownContained=false;
        }
      }
      else
      {
        notKnownContained=false;
      }
  }  // END while (notKnownContained)
  //
  // Search downwards until deepest containing volume found,
  // blocking fBlockedPhysicalVolume/BlockedReplicaNum
  //
  // 3 Cases:
  //
  // o Parameterised daughters
  //   =>Must be one G4PVParameterised daughter & voxels
  // o Positioned daughters & voxels
  // o Positioned daughters & no voxels

  noResult = true;  // noResult should be renamed to
                    // something like enteredLevel, as that is its meaning.
  do
  {
    // Determine `type' of current mother volume
    //
    targetPhysical = fHistory.GetTopVolume();
    if (!targetPhysical) { break; }
    targetLogical = targetPhysical->GetLogicalVolume();
    switch( CharacteriseDaughters(targetLogical) )
    {
      case kNormal:
        if ( targetLogical->GetVoxelHeader() )  // use optimised navigation
        {
          noResult = fvoxelNav.LevelLocate(fHistory,
                                           fBlockedPhysicalVolume,
                                           fBlockedReplicaNo,
                                           globalPoint,
                                           pGlobalDirection,
                                           considerDirection,
                                           localPoint);
        }
        else                       // do not use optimised navigation
        {
          noResult = fnormalNav.LevelLocate(fHistory,
                                            fBlockedPhysicalVolume,
                                            fBlockedReplicaNo,
                                            globalPoint,
                                            pGlobalDirection,
                                            considerDirection,
                                            localPoint);
        }
        break;
      case kReplica:
        noResult = freplicaNav.LevelLocate(fHistory,
                                           fBlockedPhysicalVolume,
                                           fBlockedReplicaNo,
                                           globalPoint,
                                           pGlobalDirection,
                                           considerDirection,
                                           localPoint);
        break;
      case kParameterised:
        if( GetDaughtersRegularStructureId(targetLogical) != 1 )
        {
          noResult = fparamNav.LevelLocate(fHistory,
                                           fBlockedPhysicalVolume,
                                           fBlockedReplicaNo,
                                           globalPoint,
                                           pGlobalDirection,
                                           considerDirection,
                                           localPoint);
        }
        else  // Regular structure
        {
          noResult = fregularNav.LevelLocate(fHistory,
                                             fBlockedPhysicalVolume,
                                             fBlockedReplicaNo,
                                             globalPoint,
                                             pGlobalDirection,
                                             considerDirection,
                                             localPoint);
        }
        break;
    }

    // LevelLocate returns true if it finds a daughter volume 
    // in which globalPoint is inside (or on the surface).

    if ( noResult )
    {
      noLevelsEntered++;
      
      // Entering a daughter after ascending
      //
      // The blocked volume is no longer valid - it was for another level
      //
      fBlockedPhysicalVolume = 0;
      fBlockedReplicaNo = -1;

      // fEntering should be false -- else blockedVolume is assumed good.
      // fEnteredDaughter is used for ExitNormal
      //
      fEntering = false;
      fEnteredDaughter = true;

      if( fExitedMother )
      {
        G4VPhysicalVolume* enteredPhysical = fHistory.GetTopVolume();
        const G4RotationMatrix* mRot = enteredPhysical->GetRotation();
        if( mRot )
        { 
          // Go deeper, i.e. move 'down' in the hierarchy
          // Apply direct rotation, not inverse
          //
          fGrandMotherExitNormal *= (*mRot);
          fChangedGrandMotherRefFrame= true;
        }
      }

#ifdef G4DEBUG_NAVIGATION
      if( fVerbose > 2 )
      { 
         G4VPhysicalVolume* enteredPhysical = fHistory.GetTopVolume();
         G4cout << "*** G4Navigator::LocateGlobalPointAndSetup() ***" << G4endl;
         G4cout << "    Entering volume: " << enteredPhysical->GetName()
                << G4endl;
      }
#endif
    }
  } while (noResult);  // Loop checking, 07.10.2016, J.Apostolakis

  fLastLocatedPointLocal = localPoint;

#ifdef G4VERBOSE
  if( fVerbose >= 4 )
  {
    G4int oldcoutPrec = G4cout.precision(8);
    G4String curPhysVol_Name("None");
    if (targetPhysical)  { curPhysVol_Name = targetPhysical->GetName(); }
    G4cout << "    Return value = new volume = " << curPhysVol_Name << G4endl;
    G4cout << "    ----- Upon exiting:" << G4endl;
    PrintState();
    if( fVerbose >= 5 )
    {
      G4cout << "Upon exiting LocateGlobalPointAndSetup():" << G4endl;
      G4cout << "    History = " << G4endl << fHistory << G4endl << G4endl;
    }
    G4cout.precision(oldcoutPrec);
  }
#endif

  fLocatedOutsideWorld= false;

  return targetPhysical;
}

// ********************************************************************
// LocateGlobalPointWithinVolume
//
// -> the state information of this Navigator and its subNavigators
//    is updated in order to start the next step at pGlobalpoint
// -> no check is performed whether pGlobalpoint is inside the 
//    original volume (this must be the case).
//
// Note: a direction could be added to the arguments, to aid in future
//       optional checking (via the old code below, flagged by OLD_LOCATE). 
//       [ This would be done only in verbose mode ]
// ********************************************************************
//
void
G4Navigator::LocateGlobalPointWithinVolume(const G4ThreeVector& pGlobalpoint)
{
#ifdef G4DEBUG_NAVIGATION
   // Check: Either step was not limited by a boundary
   //         or else the full step is no longer being taken
   assert( !fWasLimitedByGeometry );
#endif
  
   fLastLocatedPointLocal = ComputeLocalPoint(pGlobalpoint);
   fLastTriedStepComputation= false;
   fChangedGrandMotherRefFrame= false;  //  Frame for Exit Normal

#ifdef G4DEBUG_NAVIGATION
   if( fVerbose > 2 )
   { 
     G4cout << "Entering LocateGlobalWithinVolume(): History = " << G4endl;
     G4cout << fHistory << G4endl;
   }
#endif

   // For the case of Voxel (or Parameterised) volume the respective 
   // Navigator must be messaged to update its voxel information etc

   // Update the state of the Sub Navigators 
   // - in particular any voxel information they store/cache
   //
   G4VPhysicalVolume*  motherPhysical = fHistory.GetTopVolume();
   G4LogicalVolume*    motherLogical  = motherPhysical->GetLogicalVolume();
   G4SmartVoxelHeader* pVoxelHeader   = motherLogical->GetVoxelHeader();

   if ( fHistory.GetTopVolumeType()!=kReplica )
   {
     switch( CharacteriseDaughters(motherLogical) )
     {
       case kNormal:
         if ( pVoxelHeader )
         {
           fvoxelNav.VoxelLocate( pVoxelHeader, fLastLocatedPointLocal );
         }
         break;
       case kParameterised:
         if( GetDaughtersRegularStructureId(motherLogical) != 1 )
         {
           // Resets state & returns voxel node
           //
           fparamNav.ParamVoxelLocate( pVoxelHeader, fLastLocatedPointLocal );
         }
         break;
       case kReplica:
         G4Exception("G4Navigator::LocateGlobalPointWithinVolume()",
                     "GeomNav0001", FatalException,
                     "Not applicable for replicated volumes.");
         break;
     }
   }

   // Reset the state variables 
   //   - which would have been affected
   //     by the 'equivalent' call to LocateGlobalPointAndSetup
   //   - who's values have been invalidated by the 'move'.
   //
   fBlockedPhysicalVolume = 0; 
   fBlockedReplicaNo = -1;
   fEntering = false;
   fEnteredDaughter = false;  // Boundary not encountered, did not enter
   fExiting = false;
   fExitedMother = false;     // Boundary not encountered, did not exit
}

// ********************************************************************
// SetSavedState
//
// Save the state, in case this is a parasitic call
// Save fValidExitNormal, fExitNormal, fExiting, fEntering, 
//      fBlockedPhysicalVolume, fBlockedReplicaNo, fLastStepWasZero; 
// ********************************************************************
//
void G4Navigator::SetSavedState()
{
  // Note: the state of dependent objects is not currently saved.
  //   ( This means that the full state is changed by calls between
  //     SetSavedState() and RestoreSavedState(); 
  
  fSaveState.sExitNormal = fExitNormal;
  fSaveState.sValidExitNormal = fValidExitNormal;
  fSaveState.sExiting = fExiting;
  fSaveState.sEntering = fEntering;

  fSaveState.spBlockedPhysicalVolume = fBlockedPhysicalVolume;
  fSaveState.sBlockedReplicaNo = fBlockedReplicaNo, 

  fSaveState.sLastStepWasZero = fLastStepWasZero;
  
  fSaveState.sLocatedOutsideWorld = fLocatedOutsideWorld;
  fSaveState.sLastLocatedPointLocal= fLastLocatedPointLocal;
  fSaveState.sEnteredDaughter= fEnteredDaughter;
  fSaveState.sExitedMother= fExitedMother;
  fSaveState.sWasLimitedByGeometry= fWasLimitedByGeometry;

  // Even the safety sphere - if you want to change it do it explicitly!
  //
  fSaveState.sPreviousSftOrigin= fPreviousSftOrigin;
  fSaveState.sPreviousSafety= fPreviousSafety;
}

// ********************************************************************
// RestoreSavedState
//
// Restore the state (in Compute Step), in case this is a parasitic call
// ********************************************************************
//
void G4Navigator::RestoreSavedState()
{
  fExitNormal = fSaveState.sExitNormal;
  fValidExitNormal = fSaveState.sValidExitNormal;
  fExiting = fSaveState.sExiting;
  fEntering = fSaveState.sEntering;

  fBlockedPhysicalVolume = fSaveState.spBlockedPhysicalVolume;
  fBlockedReplicaNo = fSaveState.sBlockedReplicaNo, 

  fLastStepWasZero = fSaveState.sLastStepWasZero;
  
  fLocatedOutsideWorld = fSaveState.sLocatedOutsideWorld;
  fLastLocatedPointLocal= fSaveState.sLastLocatedPointLocal;
  fEnteredDaughter= fSaveState.sEnteredDaughter;
  fExitedMother= fSaveState.sExitedMother;
  fWasLimitedByGeometry= fSaveState.sWasLimitedByGeometry;
  
  // The 'expected' behaviour is to restore these too (fix 2014.05.26)
  fPreviousSftOrigin=   fSaveState.sPreviousSftOrigin;
  fPreviousSafety= fSaveState.sPreviousSafety;
}

// ********************************************************************
// ComputeStep
//
// Computes the next geometric Step: intersections with current
// mother and `daughter' volumes.
//
// NOTE:
//
// Flags on entry:
// --------------
// fValidExitNormal  - Normal of exited volume is valid (convex, not a 
//                     coincident boundary)
// fExitNormal       - Surface normal of exited volume
// fExiting          - True if have exited solid
//
// fBlockedPhysicalVolume - Ptr to exited volume (or 0)
// fBlockedReplicaNo - Replication no of exited volume
// fLastStepWasZero  - True if last Step size was almost zero.
//
// Flags on exit:
// -------------
// fValidExitNormal  - True if surface normal of exited volume is valid
// fExitNormal       - Surface normal of exited volume rotated to mothers
//                    reference system
// fExiting          - True if exiting mother
// fEntering         - True if entering `daughter' volume (or replica)
// fBlockedPhysicalVolume - Ptr to candidate (entered) volume
// fBlockedReplicaNo - Replication no of candidate (entered) volume
// fLastStepWasZero  - True if this Step size was almost zero.
// ********************************************************************
//
G4double G4Navigator::ComputeStep( const G4ThreeVector &pGlobalpoint,
                                   const G4ThreeVector &pDirection,
                                   const G4double pCurrentProposedStepLength,
                                         G4double &pNewSafety)
{
  G4ThreeVector localDirection = ComputeLocalAxis(pDirection);
  G4double Step = kInfinity;
  G4VPhysicalVolume  *motherPhysical = fHistory.GetTopVolume();
  G4LogicalVolume *motherLogical = motherPhysical->GetLogicalVolume();

  // All state relating to exiting normals must be reset
  //
  fExitNormalGlobalFrame= G4ThreeVector( 0., 0., 0.);
    // Reset value - to erase its memory
  fChangedGrandMotherRefFrame= false;
    // Reset - used for local exit normal
  fGrandMotherExitNormal= G4ThreeVector( 0., 0., 0.); 
  fCalculatedExitNormal  = false;
    // Reset for new step

  static G4ThreadLocal G4int sNavCScalls=0;
  sNavCScalls++;

  fLastTriedStepComputation= true; 

#ifdef G4VERBOSE
  if( fVerbose > 0 )
  {
    G4cout << "*** G4Navigator::ComputeStep: ***" << G4endl; 
    G4cout << "    Volume = " << motherPhysical->GetName() 
           << " - Proposed step length = " << pCurrentProposedStepLength
           << G4endl; 
#ifdef G4DEBUG_NAVIGATION
    if( fVerbose >= 2 )
    {
      G4cout << "  Called with the arguments: " << G4endl
             << "  Globalpoint = " << std::setw(25) << pGlobalpoint << G4endl
             << "  Direction   = " << std::setw(25) << pDirection << G4endl;
      if( fVerbose >= 4 )
      {
        G4cout << "  ---- Upon entering : State" << G4endl;
        PrintState();
      }
    }
#endif
  }
#endif

  G4ThreeVector newLocalPoint = ComputeLocalPoint(pGlobalpoint);
  if( newLocalPoint != fLastLocatedPointLocal )
  {
    // Check whether the relocation is within safety
    //
    G4ThreeVector oldLocalPoint = fLastLocatedPointLocal;
    G4double moveLenSq = (newLocalPoint-oldLocalPoint).mag2();

    if ( moveLenSq >= fSqTol )
    {
#ifdef G4VERBOSE
      ComputeStepLog(pGlobalpoint, moveLenSq);
#endif
      // Relocate the point within the same volume
      //
      LocateGlobalPointWithinVolume( pGlobalpoint );
      fLastTriedStepComputation= true;     // Ensure that this is set again !!
    }
  }
  if ( fHistory.GetTopVolumeType()!=kReplica )
  {
    switch( CharacteriseDaughters(motherLogical) )
    {
      case kNormal:
        if ( motherLogical->GetVoxelHeader() )
        {
          Step = fvoxelNav.ComputeStep(fLastLocatedPointLocal,
                                       localDirection,
                                       pCurrentProposedStepLength,
                                       pNewSafety,
                                       fHistory,
                                       fValidExitNormal,
                                       fExitNormal,
                                       fExiting,
                                       fEntering,
                                       &fBlockedPhysicalVolume,
                                       fBlockedReplicaNo);
      
        }
        else
        {
          if( motherPhysical->GetRegularStructureId() == 0 )
          {
            Step = fnormalNav.ComputeStep(fLastLocatedPointLocal,
                                          localDirection,
                                          pCurrentProposedStepLength,
                                          pNewSafety,
                                          fHistory,
                                          fValidExitNormal,
                                          fExitNormal,
                                          fExiting,
                                          fEntering,
                                          &fBlockedPhysicalVolume,
                                          fBlockedReplicaNo);
          }
          else  // Regular (non-voxelised) structure
          {
            LocateGlobalPointAndSetup( pGlobalpoint, &pDirection, true, true );
            fLastTriedStepComputation= true; // Ensure that this is set again !!
            //
            // if physical process limits the step, the voxel will not be the
            // one given by ComputeStepSkippingEqualMaterials() and the local
            // point will be wrongly calculated.

            // There is a problem: when msc limits the step and the point is
            // assigned wrongly to phantom in previous step (while it is out
            // of the container volume). Then LocateGlobalPointAndSetup() has
            // reset the history topvolume to world.
            //
            if(fHistory.GetTopVolume()->GetRegularStructureId() == 0 )
            { 
              G4Exception("G4Navigator::ComputeStep()",
                          "GeomNav1001", JustWarning,
                "Point is relocated in voxels, while it should be outside!");
              Step = fnormalNav.ComputeStep(fLastLocatedPointLocal,
                                            localDirection,
                                            pCurrentProposedStepLength,
                                            pNewSafety,
                                            fHistory,
                                            fValidExitNormal,
                                            fExitNormal,
                                            fExiting,
                                            fEntering,
                                            &fBlockedPhysicalVolume,
                                            fBlockedReplicaNo);
            }
            else
            {
              Step = fregularNav.
                   ComputeStepSkippingEqualMaterials(fLastLocatedPointLocal,
                                                     localDirection,
                                                     pCurrentProposedStepLength,
                                                     pNewSafety,
                                                     fHistory,
                                                     fValidExitNormal,
                                                     fExitNormal,
                                                     fExiting,
                                                     fEntering,
                                                     &fBlockedPhysicalVolume,
                                                     fBlockedReplicaNo,
                                                     motherPhysical);
            }
          }
        }
        break;
      case kParameterised:
        if( GetDaughtersRegularStructureId(motherLogical) != 1 )
        {
          Step = fparamNav.ComputeStep(fLastLocatedPointLocal,
                                       localDirection,
                                       pCurrentProposedStepLength,
                                       pNewSafety,
                                       fHistory,
                                       fValidExitNormal,
                                       fExitNormal,
                                       fExiting,
                                       fEntering,
                                       &fBlockedPhysicalVolume,
                                       fBlockedReplicaNo);
        }
        else  // Regular structure
        {
          Step = fregularNav.ComputeStep(fLastLocatedPointLocal,
                                         localDirection,
                                         pCurrentProposedStepLength,
                                         pNewSafety,
                                         fHistory,
                                         fValidExitNormal,
                                         fExitNormal,
                                         fExiting,
                                         fEntering,
                                         &fBlockedPhysicalVolume,
                                         fBlockedReplicaNo);
        }
        break;
      case kReplica:
        G4Exception("G4Navigator::ComputeStep()", "GeomNav0001",
                    FatalException, "Not applicable for replicated volumes.");
        break;
    }
  }
  else
  {
    // In the case of a replica, it must handle the exiting
    // edge/corner problem by itself
    //
    G4bool exitingReplica = fExitedMother;
    G4bool calculatedExitNormal;
    Step = freplicaNav.ComputeStep(pGlobalpoint,
                                   pDirection,
                                   fLastLocatedPointLocal,
                                   localDirection,
                                   pCurrentProposedStepLength,
                                   pNewSafety,
                                   fHistory,
                                   fValidExitNormal,
                                   calculatedExitNormal,
                                   fExitNormal,
                                   exitingReplica,
                                   fEntering,
                                   &fBlockedPhysicalVolume,
                                   fBlockedReplicaNo);
    fExiting= exitingReplica;                          
    fCalculatedExitNormal= calculatedExitNormal;
  }

  // Remember last safety origin & value.
  //
  fPreviousSftOrigin = pGlobalpoint;
  fPreviousSafety = pNewSafety; 

  // Count zero steps - one can occur due to changing momentum at a boundary
  //                  - one, two (or a few) can occur at common edges between
  //                    volumes
  //                  - more than two is likely a problem in the geometry
  //                    description or the Navigation 

  // Rule of thumb: likely at an Edge if two consecutive steps are zero,
  //                because at least two candidate volumes must have been
  //                checked
  //
  fLocatedOnEdge   = fLastStepWasZero && (Step==0.0);
  fLastStepWasZero = (Step<fMinStep);
  if (fPushed)  { fPushed = fLastStepWasZero; }

  // Handle large number of consecutive zero steps
  //
  if ( fLastStepWasZero )
  {
    fNumberZeroSteps++;
#ifdef G4DEBUG_NAVIGATION
    if( fNumberZeroSteps > 1 )
    {
       G4cout << "G4Navigator::ComputeStep(): another 'zero' step, # "
              << fNumberZeroSteps
              << ", at " << pGlobalpoint
              << ", in volume " << motherPhysical->GetName()
              << ", nav-comp-step calls # " << sNavCScalls
              << ", Step= " << Step
              << G4endl;
    }
#endif
    if( fNumberZeroSteps > fActionThreshold_NoZeroSteps-1 )
    {
       // Act to recover this stuck track. Pushing it along direction
       //
       Step += 100*kCarTolerance;
#ifdef G4VERBOSE
       if ((!fPushed) && (fWarnPush))
       {
         std::ostringstream message;
         message.precision(16);
         message << "Track stuck or not moving." << G4endl
                 << "          Track stuck, not moving for " 
                 << fNumberZeroSteps << " steps" << G4endl
                 << "          in volume -" << motherPhysical->GetName()
                 << "- at point " << pGlobalpoint
                 << " (local point " << newLocalPoint << ")" << G4endl
                 << "          direction: " << pDirection
                 << " (local direction: " << localDirection << ")." << G4endl
                 << "          Potential geometry or navigation problem !"
                 << G4endl
                 << "          Trying pushing it of " << Step << " mm ...";
         G4Exception("G4Navigator::ComputeStep()", "GeomNav1002",
                     JustWarning, message, "Potential overlap in geometry!");
       }
#endif
       fPushed = true;
    }
    if( fNumberZeroSteps > fAbandonThreshold_NoZeroSteps-1 )
    {
      // Must kill this stuck track
      //
      std::ostringstream message;
      message << "Stuck Track: potential geometry or navigation problem."
              << G4endl
              << "        Track stuck, not moving for " 
              << fNumberZeroSteps << " steps" << G4endl
              << "        in volume -" << motherPhysical->GetName()
              << "- at point " << pGlobalpoint << G4endl
              << "        direction: " << pDirection << ".";
#ifdef G4VERBOSE
      if ( fWarnPush )
      {
         motherPhysical->CheckOverlaps(5000, false);
      }
#endif
      G4Exception("G4Navigator::ComputeStep()", "GeomNav0003",
                  EventMustBeAborted, message);
    }
  }
  else
  {
    if (!fPushed)  fNumberZeroSteps = 0;
  }

  fEnteredDaughter = fEntering;   // I expect to enter a volume in this Step
  fExitedMother = fExiting;

  fStepEndPoint = pGlobalpoint
                + std::min(Step,pCurrentProposedStepLength) * pDirection;
  fLastStepEndPointLocal = fLastLocatedPointLocal + Step * localDirection; 

  if( fExiting )
  {
#ifdef G4DEBUG_NAVIGATION
    if( fVerbose > 2 )
    { 
      G4cout << " At G4Nav CompStep End - if(exiting) - fExiting= " << fExiting 
             << " fValidExitNormal = " << fValidExitNormal  << G4endl;
      G4cout << " fExitNormal= " << fExitNormal << G4endl;
    }
#endif

    if(fValidExitNormal || fCalculatedExitNormal)
    {
      if (  fHistory.GetTopVolumeType()!=kReplica )
      {
        // Convention: fExitNormal is in the 'grand-mother' coordinate system
        //
        fGrandMotherExitNormal= fExitNormal;
        fCalculatedExitNormal= true;
      }
      else
      {
        fGrandMotherExitNormal = fExitNormal;
      }
    }
    else
    {  
      // We must calculate the normal anyway (in order to have it if requested)
      //
      G4ThreeVector finalLocalPoint =
            fLastLocatedPointLocal + localDirection*Step;

      if (  fHistory.GetTopVolumeType()!=kReplica )
      {
        // Find normal in the 'mother' coordinate system
        //
        G4ThreeVector exitNormalMotherFrame=
           motherLogical->GetSolid()->SurfaceNormal(finalLocalPoint);
        
        // Transform it to the 'grand-mother' coordinate system
        //
        const G4RotationMatrix* mRot = motherPhysical->GetRotation();
        if( mRot )
        {
          fChangedGrandMotherRefFrame= true;           
          fGrandMotherExitNormal = (*mRot).inverse() * exitNormalMotherFrame;
        }
        else
        {
          fGrandMotherExitNormal = exitNormalMotherFrame;
        }

        // Do not set fValidExitNormal -- this signifies
        // that the solid is convex!
        //
        fCalculatedExitNormal= true;
      }
      else
      {
        fCalculatedExitNormal = false;
        //
        // Nothing can be done at this stage currently - to solve this
        // Replica Navigation must have calculated the normal for this case
        // already.
        // Cases: mother is not convex, and exit is at previous replica level

#ifdef G4DEBUG_NAVIGATION
        G4ExceptionDescription desc;

        desc << "Problem in ComputeStep:  Replica Navigation did not provide"
             << " valid exit Normal. " << G4endl;
        desc << " Do not know how calculate it in this case." << G4endl;
        desc << "  Location    = " << finalLocalPoint << G4endl;
        desc << "  Volume name = " << motherPhysical->GetName()
             << "  copy/replica No = " << motherPhysical->GetCopyNo() << G4endl;
        G4Exception("G4Navigator::ComputeStep()", "GeomNav0003",
                    JustWarning, desc, "Normal not available for exiting.");
#endif
      }
    }

    // Now transform it to the global reference frame !!
    //
    if( fValidExitNormal || fCalculatedExitNormal )
    {
      G4int depth= fHistory.GetDepth();
      if( depth > 0 )
      {
        G4AffineTransform GrandMotherToGlobalTransf =
          fHistory.GetTransform(depth-1).Inverse();
        fExitNormalGlobalFrame =
          GrandMotherToGlobalTransf.TransformAxis( fGrandMotherExitNormal );
      }
      else
      {
        fExitNormalGlobalFrame= fGrandMotherExitNormal;
      }
    }
    else
    {
      fExitNormalGlobalFrame= G4ThreeVector( 0., 0., 0.);
    }
  }

  if( (Step == pCurrentProposedStepLength) && (!fExiting) && (!fEntering) )
  {
    // This if Step is not really limited by the geometry.
    // The Navigator is obliged to return "infinity"
    //
    Step = kInfinity;
  }

#ifdef G4VERBOSE
  if( fVerbose > 1 )
  {
    if( fVerbose >= 4 )
    {
      G4cout << "    ----- Upon exiting :" << G4endl;
      PrintState();
    }
    G4cout << "  Returned step= " << Step;
    if( fVerbose > 5 )   G4cout << G4endl;
    if( Step == kInfinity )
    {
       G4cout << " Requested step= " << pCurrentProposedStepLength ;
       if( fVerbose > 5) G4cout << G4endl;
    }
    G4cout << "  Safety = " << pNewSafety << G4endl;
  }
#endif

  return Step;
}

// ********************************************************************
// CheckNextStep
//
// Compute the step without altering the navigator state
// ********************************************************************
//
G4double G4Navigator::CheckNextStep( const G4ThreeVector& pGlobalpoint,
                                     const G4ThreeVector& pDirection,
                                     const G4double pCurrentProposedStepLength,
                                           G4double& pNewSafety)
{
  G4double step;

  // Save the state, for this parasitic call
  //
  SetSavedState();

  step = ComputeStep ( pGlobalpoint, 
                       pDirection,
                       pCurrentProposedStepLength, 
                       pNewSafety ); 

  // It is a parasitic call, so attempt to restore the key parts of the state
  //
  RestoreSavedState(); 
  // NOTE: the state of the current subnavigator is NOT restored.
  // ***> TODO: restore subnavigator state
  //            if( last_located)       Need Position of last location
  //            if( last_computed step) Need Endposition of last step
  
  return step; 
}

// ********************************************************************
// ResetState
//
// Resets stack and minimum of navigator state `machine'
// ********************************************************************
//
void G4Navigator::ResetState()
{
  fWasLimitedByGeometry  = false;
  fEntering              = false;
  fExiting               = false;
  fLocatedOnEdge         = false;
  fLastStepWasZero       = false;
  fEnteredDaughter       = false;
  fExitedMother          = false;
  fPushed                = false;

  fValidExitNormal       = false;
  fChangedGrandMotherRefFrame= false;
  fCalculatedExitNormal  = false;

  fExitNormal            = G4ThreeVector(0,0,0);
  fGrandMotherExitNormal = G4ThreeVector(0,0,0);
  fExitNormalGlobalFrame = G4ThreeVector(0,0,0);

  fPreviousSftOrigin     = G4ThreeVector(0,0,0);
  fPreviousSafety        = 0.0; 

  fNumberZeroSteps       = 0;
    
  fBlockedPhysicalVolume = 0;
  fBlockedReplicaNo      = -1;

  fLastLocatedPointLocal = G4ThreeVector( kInfinity, -kInfinity, 0.0 ); 
  fLocatedOutsideWorld   = false;
}

// ********************************************************************
// SetupHierarchy
//
// Renavigates & resets hierarchy described by current history
// o Reset volumes
// o Recompute transforms and/or solids of replicated/parameterised volumes
// ********************************************************************
//
void G4Navigator::SetupHierarchy()
{
  G4int i;
  const G4int cdepth = fHistory.GetDepth();
  G4VPhysicalVolume *current;
  G4VSolid *pSolid;
  G4VPVParameterisation *pParam;

  for ( i=1; i<=cdepth; i++ )
  {
    current = fHistory.GetVolume(i);
    switch ( fHistory.GetVolumeType(i) )
    {
      case kNormal:
        break;
      case kReplica:
        freplicaNav.ComputeTransformation(fHistory.GetReplicaNo(i), current);
        break;
      case kParameterised:
        G4int replicaNo;
        pParam = current->GetParameterisation();
        replicaNo = fHistory.GetReplicaNo(i);
        pSolid = pParam->ComputeSolid(replicaNo, current);

        // Set up dimensions & transform in solid/physical volume
        //
        pSolid->ComputeDimensions(pParam, replicaNo, current);
        pParam->ComputeTransformation(replicaNo, current);

        G4TouchableHistory *pTouchable= 0;
        if( pParam->IsNested() )
        {
           pTouchable= new G4TouchableHistory( fHistory );
           pTouchable->MoveUpHistory(); // Move up to the parent level 
             // Adequate only if Nested at the Branch level (last)
           // To extend to other cases:  
           // pTouchable->MoveUpHistory(cdepth-i-1);
             // Move to the parent level of *Current* level
             // Could replace this line and constructor with a revised
             // c-tor for History(levels to drop)
        }
        // Set up the correct solid and material in Logical Volume
        //
        G4LogicalVolume *pLogical = current->GetLogicalVolume();
        pLogical->SetSolid( pSolid );
        pLogical->UpdateMaterial( pParam ->
          ComputeMaterial(replicaNo, current, pTouchable) );
        delete pTouchable;
        break;
    }
  }
}

// ********************************************************************
// GetLocalExitNormal
//
// Obtains the Normal vector to a surface (in local coordinates)
// pointing out of previous volume and into current volume
// ********************************************************************
//
G4ThreeVector G4Navigator::GetLocalExitNormal( G4bool* valid )
{
  G4ThreeVector    ExitNormal(0.,0.,0.);
  G4VSolid        *currentSolid=0;
  G4LogicalVolume *candidateLogical;
  if ( fLastTriedStepComputation ) 
  {
    // use fLastLocatedPointLocal and next candidate volume
    //
    G4ThreeVector nextSolidExitNormal(0.,0.,0.);

    if( fEntering && (fBlockedPhysicalVolume!=0) ) 
    { 
      candidateLogical= fBlockedPhysicalVolume->GetLogicalVolume();
      if( candidateLogical ) 
      {
        // fLastStepEndPointLocal is in the coordinates of the mother
        // we need it in the daughter's coordinate system.

        // The following code should also work in case of Replica
        {
          // First transform fLastLocatedPointLocal to the new daughter
          // coordinates
          //
          G4AffineTransform MotherToDaughterTransform=
            GetMotherToDaughterTransform( fBlockedPhysicalVolume, 
                                          fBlockedReplicaNo,
                                          VolumeType(fBlockedPhysicalVolume) ); 
          G4ThreeVector daughterPointOwnLocal= 
            MotherToDaughterTransform.TransformPoint( fLastStepEndPointLocal ); 

          // OK if it is a parameterised volume
          //
          EInside  inSideIt; 
          G4bool   onSurface;
          G4double safety= -1.0; 
          currentSolid= candidateLogical->GetSolid(); 
          inSideIt  =  currentSolid->Inside(daughterPointOwnLocal); 
          onSurface =  (inSideIt == kSurface); 
          if( ! onSurface ) 
          {
            if( inSideIt == kOutside )
            {
              safety = (currentSolid->DistanceToIn(daughterPointOwnLocal)); 
              onSurface = safety < 100.0 * kCarTolerance; 
            }
            else if (inSideIt == kInside )
            {
              safety = (currentSolid->DistanceToOut(daughterPointOwnLocal)); 
              onSurface = safety < 100.0 * kCarTolerance; 
            }
          }

          if( onSurface ) 
          {
            nextSolidExitNormal =
              currentSolid->SurfaceNormal(daughterPointOwnLocal); 
 
            // Entering the solid ==> opposite
            //
            ExitNormal = -nextSolidExitNormal;
            fCalculatedExitNormal= true;
          }
          else
          {
#ifdef G4VERBOSE
            if(( fVerbose == 1 ) && ( fCheck ))
            {
              std::ostringstream message;
              message << "Point not on surface ! " << G4endl
                      << "  Point           = "
                      << daughterPointOwnLocal << G4endl 
                      << "  Physical volume = "
                      << fBlockedPhysicalVolume->GetName() << G4endl
                      << "  Logical volume  = "
                      << candidateLogical->GetName() << G4endl
                      << "  Solid           = " << currentSolid->GetName() 
                      << "  Type            = "
                      << currentSolid->GetEntityType() << G4endl
                      << *currentSolid << G4endl;
              if( inSideIt == kOutside )
              { 
                message << "Point is Outside. " << G4endl
                        << "  Safety (from outside) = " << safety << G4endl;
              }
              else // if( inSideIt == kInside ) 
              {
                message << "Point is Inside. " << G4endl
                        << "  Safety (from inside) = " << safety << G4endl;
              }
              G4Exception("G4Navigator::GetLocalExitNormal()", "GeomNav1001",
                          JustWarning, message);
            }
#endif
          }
          *valid = onSurface;   //   was =true;
        }
      }
    }
    else if ( fExiting ) 
    {
      ExitNormal = fGrandMotherExitNormal;
      *valid = true;
      fCalculatedExitNormal= true;  // Should be true already
    }
    else  // i.e.  ( fBlockedPhysicalVolume == 0 )
    {
      *valid = false;
      G4Exception("G4Navigator::GetLocalExitNormal()",
                  "GeomNav0003", JustWarning, 
                  "Incorrect call to GetLocalSurfaceNormal." );
    }
  }
  else //  ( ! fLastTriedStepComputation ) ie. last call was to Locate
  {
    if ( EnteredDaughterVolume() )
    {
      G4VSolid* daughterSolid =fHistory.GetTopVolume()->GetLogicalVolume()
                                                      ->GetSolid();
      ExitNormal= -(daughterSolid->SurfaceNormal(fLastLocatedPointLocal));
      if( std::fabs(ExitNormal.mag2()-1.0 ) > CLHEP::perMillion )
      {
        G4ExceptionDescription desc;
        desc << " Parameters of solid: " << *daughterSolid
             << " Point for surface = " << fLastLocatedPointLocal << std::endl;
        G4Exception("G4Navigator::GetLocalExitNormal()",
                    "GeomNav0003", FatalException, desc,
                    "Surface Normal returned by Solid is not a Unit Vector." );
      }
      fCalculatedExitNormal= true;
      *valid = true;
    }
    else
    {
      if( fExitedMother )
      {
        ExitNormal = fGrandMotherExitNormal;
        *valid = true;
        fCalculatedExitNormal= true;
      }
      else  // We are not at a boundary. ExitNormal remains (0,0,0)
      { 
        *valid = false;
        fCalculatedExitNormal= false; 
        G4ExceptionDescription message; 
        message << "Function called when *NOT* at a Boundary." << G4endl;
        message << "Exit Normal not calculated." << G4endl;
        G4Exception("G4Navigator::GetLocalExitNormal()",
                    "GeomNav0003", JustWarning, message); 
      }
    }
  }
  return ExitNormal;
}

// ********************************************************************
// GetMotherToDaughterTransform
//
// Obtains the mother to daughter affine transformation
// ********************************************************************
//
G4AffineTransform
G4Navigator::GetMotherToDaughterTransform( G4VPhysicalVolume *pEnteringPhysVol,
                                           G4int   enteringReplicaNo,
                                           EVolume enteringVolumeType ) 
{
  switch (enteringVolumeType)
  {
    case kNormal:  // Nothing is needed to prepare the transformation
      break;       // It is stored already in the physical volume (placement)
    case kReplica: // Sets the transform in the Replica - tbc
      G4Exception("G4Navigator::GetMotherToDaughterTransform()",
                  "GeomNav0001", FatalException,
                  "Method NOT Implemented yet for replica volumes.");
      break;
    case kParameterised:
      if( pEnteringPhysVol->GetRegularStructureId() == 0 )
      {
        G4VPVParameterisation *pParam =
          pEnteringPhysVol->GetParameterisation();
        G4VSolid* pSolid =
          pParam->ComputeSolid(enteringReplicaNo, pEnteringPhysVol);
        pSolid->ComputeDimensions(pParam, enteringReplicaNo, pEnteringPhysVol);

        // Sets the transform in the Parameterisation
        //
        pParam->ComputeTransformation(enteringReplicaNo, pEnteringPhysVol);

        // Set the correct solid and material in Logical Volume
        //
        G4LogicalVolume* pLogical = pEnteringPhysVol->GetLogicalVolume();
        pLogical->SetSolid( pSolid );
      }
      break;
  }
  return G4AffineTransform(pEnteringPhysVol->GetRotation(), 
                           pEnteringPhysVol->GetTranslation()).Invert(); 
}


// ********************************************************************
// GetLocalExitNormalAndCheck
//
// Obtains the Normal vector to a surface (in local coordinates)
// pointing out of previous volume and into current volume, and
// checks the current point against expected 'local' value.
// ********************************************************************
//
G4ThreeVector
G4Navigator::GetLocalExitNormalAndCheck( 
#ifdef G4DEBUG_NAVIGATION
                           const G4ThreeVector& ExpectedBoundaryPointGlobal,
#else
                           const G4ThreeVector&,
#endif
                                 G4bool*        pValid)
{
#ifdef G4DEBUG_NAVIGATION
  // Check Current point against expected 'local' value
  //
  if ( fLastTriedStepComputation ) 
  {
    G4ThreeVector ExpectedBoundaryPointLocal;

    const G4AffineTransform& GlobalToLocal= GetGlobalToLocalTransform(); 
    ExpectedBoundaryPointLocal =
      GlobalToLocal.TransformPoint( ExpectedBoundaryPointGlobal ); 

    // Add here:  Comparison against expected position,
    //            i.e. the endpoint of ComputeStep
  }
#endif
  
  return GetLocalExitNormal( pValid); 
}


// ********************************************************************
// GetGlobalExitNormal
//
// Obtains the Normal vector to a surface (in global coordinates)
// pointing out of previous volume and into current volume
// ********************************************************************
//
G4ThreeVector 
G4Navigator::GetGlobalExitNormal(const G4ThreeVector& IntersectPointGlobal,
                                       G4bool*        pNormalCalculated)
{
  G4bool         validNormal;
  G4ThreeVector  localNormal, globalNormal;
 
  G4bool usingStored = fCalculatedExitNormal && (
       ( fLastTriedStepComputation && fExiting ) // Just calculated it
       ||                                        // No locate in between
       ( !fLastTriedStepComputation
          && (IntersectPointGlobal-fStepEndPoint).mag2() < 10.0*fSqTol ) );
           // Calculated it 'just' before & then called locate
           // but it did not move position
  
  if( usingStored )
  {
    // This was computed in last call to ComputeStep
    // and only if it arrived at boundary
    //
    globalNormal = fExitNormalGlobalFrame; 
    G4double  normMag2 = globalNormal.mag2(); 
    if( std::fabs ( normMag2 - 1.0 ) < perMillion )  // Value is good 
    {
       *pNormalCalculated = true; // ComputeStep always computes it if Exiting
                                  // (fExiting==true)
    }
    else
    {
       G4ExceptionDescription message;
       message.precision(10); 
       message << " WARNING> Expected normal-global-frame to be valid, "
               << " i.e. a unit vector!" << G4endl
               << "  - but |normal|   = "  << std::sqrt(normMag2)
               << "  - and |normal|^2 = "  << normMag2 << G4endl
               << " which differs from 1.0 by " << normMag2 - 1.0 << G4endl
               << "   n = " << fExitNormalGlobalFrame << G4endl
               << " Global point: " << IntersectPointGlobal << G4endl
               << " Volume: " << fHistory.GetTopVolume()->GetName() << G4endl;
#ifdef G4VERBOSE
       G4LogicalVolume* candLog = fHistory.GetTopVolume()->GetLogicalVolume();
       if ( candLog )
       {
         message << " Solid: " << candLog->GetSolid()->GetName()
                 << ", Type: " << candLog->GetSolid()->GetEntityType() << G4endl
                 << *candLog->GetSolid() << G4endl;
       }
#endif
       message << "============================================================"
               << G4endl;
       G4int oldVerbose = fVerbose; 
       fVerbose=4; 
       message << "   State of Navigator: " << G4endl;
       message << *this << G4endl;
       fVerbose = oldVerbose; 
       message << "============================================================"
               << G4endl;

       G4Exception("G4Navigator::GetGlobalExitNormal()",
                   "GeomNav0003",JustWarning, message,
              "Value obtained from stored global-normal is not a unit vector.");

       // (Re)Compute it now -- as either it was not computed, or it is wrong.
       //
       localNormal = GetLocalExitNormalAndCheck(IntersectPointGlobal,
                                                &validNormal);
       *pNormalCalculated = fCalculatedExitNormal;

       G4AffineTransform localToGlobal = GetLocalToGlobalTransform();
       globalNormal = localToGlobal.TransformAxis( localNormal );
    }
  }
  else
  {
    localNormal = GetLocalExitNormalAndCheck(IntersectPointGlobal,&validNormal);
    *pNormalCalculated = fCalculatedExitNormal;

#ifdef G4DEBUG_NAVIGATION
    usingStored= false;

    if( (!validNormal) && !fCalculatedExitNormal )
    {
      G4ExceptionDescription edN;
      edN << "  Calculated = " << fCalculatedExitNormal << G4endl;
      edN << "   Entering= "  << fEntering << G4endl;
      G4int oldVerbose= this->GetVerboseLevel();
      this->SetVerboseLevel(4);
      edN << "   State of Navigator: " << G4endl;
      edN << *this << G4endl;
      this->SetVerboseLevel( oldVerbose );
       
      G4Exception("G4Navigator::GetGlobalExitNormal()",
                  "GeomNav0003", JustWarning, edN,
                  "LocalExitNormalAndCheck() did not calculate Normal.");
     }
#endif
     
     G4double localMag2= localNormal.mag2();
     if( validNormal && (std::fabs(localMag2-1.0)) > CLHEP::perMillion )
     {
       G4ExceptionDescription edN;
       edN.precision(10); 
       edN << "G4Navigator::GetGlobalExitNormal: "
           << "  Using Local Normal - from call to GetLocalExitNormalAndCheck. "
           << G4endl
           << "  Local  Exit Normal : " << " || = " << std::sqrt(localMag2) 
           << " vec = " << localNormal << G4endl
           << "  Global Exit Normal : " << " || = " << globalNormal.mag() 
           << " vec = " << globalNormal << G4endl
           << "  Global point: " << IntersectPointGlobal << G4endl;
       edN << "  Calculated It      = " << fCalculatedExitNormal << G4endl
           << "  Volume: " << fHistory.GetTopVolume()->GetName() << G4endl;
#ifdef G4VERBOSE
       G4LogicalVolume* candLog = fHistory.GetTopVolume()->GetLogicalVolume();
       if ( candLog )
       {
         edN << "  Solid: " << candLog->GetSolid()->GetName()
             << ", Type: " << candLog->GetSolid()->GetEntityType() << G4endl
             << *candLog->GetSolid();
       }
#endif
       G4Exception("G4Navigator::GetGlobalExitNormal()",
                   "GeomNav0003",JustWarning, edN,
                   "Value obtained from new local *solid* is incorrect.");
       localNormal = localNormal.unit(); // Should we correct it ??
     }
     G4AffineTransform localToGlobal = GetLocalToGlobalTransform();
     globalNormal = localToGlobal.TransformAxis( localNormal );
  }

#ifdef G4DEBUG_NAVIGATION
  if( usingStored )
  {
    G4ThreeVector globalNormAgn; 

    localNormal= GetLocalExitNormalAndCheck(IntersectPointGlobal, &validNormal);
    
    G4AffineTransform localToGlobal = GetLocalToGlobalTransform();
    globalNormAgn = localToGlobal.TransformAxis( localNormal );
    
    // Check the value computed against fExitNormalGlobalFrame
    G4ThreeVector diffNorm = globalNormAgn - fExitNormalGlobalFrame;
    if( diffNorm.mag2() > perMillion*CLHEP::perMillion)
    {
      G4ExceptionDescription edDfn;
      edDfn << "Found difference in normals in case of exiting mother "
            << "- when Get is called after ComputingStep " << G4endl;
      edDfn << "  Magnitude of diff =      " << diffNorm.mag() << G4endl;
      edDfn << "  Normal stored (Global)     = " << fExitNormalGlobalFrame
            << G4endl;
      edDfn << "  Global Computed from Local = " << globalNormAgn << G4endl;
      G4Exception("G4Navigator::GetGlobalExitNormal()", "GeomNav0003",
                  JustWarning, edDfn);
    }
  }
#endif

  // Synchronise stored global exit normal as possibly re-computed here
  //
  fExitNormalGlobalFrame = globalNormal;

  return globalNormal;
}

// To make the new Voxel Safety the default, uncomment the next line
#define  G4NEW_SAFETY  1

// ********************************************************************
// ComputeSafety
//
// It assumes that it will be 
//  i) called at the Point in the same volume as the EndPoint of the
//     ComputeStep.
// ii) after (or at the end of) ComputeStep OR after the relocation.
// ********************************************************************
//
G4double G4Navigator::ComputeSafety( const G4ThreeVector &pGlobalpoint,
                                     const G4double pMaxLength,
                                     const G4bool keepState)
{
  G4double newSafety = 0.0;

#ifdef G4DEBUG_NAVIGATION
  G4int oldcoutPrec = G4cout.precision(8);
  if( fVerbose > 0 )
  {
    G4cout << "*** G4Navigator::ComputeSafety: ***" << G4endl
           << "    Called at point: " << pGlobalpoint << G4endl;

    G4VPhysicalVolume  *motherPhysical = fHistory.GetTopVolume();
    G4cout << "    Volume = " << motherPhysical->GetName() 
           << " - Maximum length = " << pMaxLength << G4endl; 
    if( fVerbose >= 4 )
    {
       G4cout << "    ----- Upon entering Compute Safety:" << G4endl;
       PrintState();
    }
  }
#endif

  G4double distEndpointSq = (pGlobalpoint-fStepEndPoint).mag2(); 
  G4bool   stayedOnEndpoint  = distEndpointSq < kCarTolerance*kCarTolerance; 
  G4bool   endpointOnSurface = fEnteredDaughter || fExitedMother;

  if( endpointOnSurface && stayedOnEndpoint )
    {
#ifdef G4DEBUG_NAVIGATION
      if( fVerbose >= 2 )
      {
        G4cout << "    G4Navigator::ComputeSafety() finds that point - "
        << pGlobalpoint << " - is on surface " << G4endl;
        if( fEnteredDaughter ) { G4cout << "   entered new daughter volume"; }
        if( fExitedMother )    { G4cout << "   and exited previous volume."; }
        G4cout << G4endl;
        G4cout << " EndPoint was = " << fStepEndPoint << G4endl;
      }
#endif
      newSafety = 0.0;
      // return newSafety;
    }
  else // if( !(endpointOnSurface && stayedOnEndpoint) )
  {
    if (keepState)  { SetSavedState(); }
    
    // Pseudo-relocate to this point (updates voxel information only)
    //
    LocateGlobalPointWithinVolume( pGlobalpoint ); 
      // --->> DANGER: Side effects on sub-navigator voxel information <<---
      //       Could be replaced again by 'granular' calls to sub-navigator
      //       locates (similar side-effects, but faster.  
      //       Solutions:
      //        1) Re-locate (to where?)
      //        2) Insure that the methods using (G4ComputeStep?)
      //           does a relocation (if information is disturbed only ?)

#ifdef G4DEBUG_NAVIGATION
    if( fVerbose >= 2 )
    {
      G4cout << "  G4Navigator::ComputeSafety() relocates-in-volume to point: "
             << pGlobalpoint << G4endl;
    }
#endif 
    G4VPhysicalVolume *motherPhysical = fHistory.GetTopVolume();
    G4LogicalVolume *motherLogical = motherPhysical->GetLogicalVolume();
    G4SmartVoxelHeader* pVoxelHeader = motherLogical->GetVoxelHeader();
    G4ThreeVector localPoint = ComputeLocalPoint(pGlobalpoint);

    if ( fHistory.GetTopVolumeType()!=kReplica )
    {
      switch(CharacteriseDaughters(motherLogical))
      {
        case kNormal:
          if ( pVoxelHeader )
          {
#ifdef G4NEW_SAFETY
            G4double safetyTwo = fpVoxelSafety->ComputeSafety(localPoint,
                                           *motherPhysical, pMaxLength);
            newSafety= safetyTwo;   // Faster and best available
#else
            G4double safetyOldVoxel;
            safetyOldVoxel =
              fvoxelNav.ComputeSafety(localPoint,fHistory,pMaxLength);
            newSafety= safetyOldVoxel;
#endif
          }
          else
          {
            newSafety=fnormalNav.ComputeSafety(localPoint,fHistory,pMaxLength);
          }
          break;
        case kParameterised:
          if( GetDaughtersRegularStructureId(motherLogical) != 1 )
          {
            newSafety=fparamNav.ComputeSafety(localPoint,fHistory,pMaxLength);
          }
          else  // Regular structure
          {
            newSafety=fregularNav.ComputeSafety(localPoint,fHistory,pMaxLength);
          }
          break;
        case kReplica:
          G4Exception("G4Navigator::ComputeSafety()", "GeomNav0001",
                      FatalException, "Not applicable for replicated volumes.");
          break;
      }
    }
    else
    {
      newSafety = freplicaNav.ComputeSafety(pGlobalpoint, localPoint,
                                            fHistory, pMaxLength);
    }
    
    if (keepState)
    {
      RestoreSavedState();
      // This now overwrites the values of the Safety 'sphere' (correction)
    }
    
    // Remember last safety origin & value
    //
    // We overwrite the Safety 'sphere' - keeping old behaviour
    fPreviousSftOrigin = pGlobalpoint;
    fPreviousSafety = newSafety;
  }
  
#ifdef G4DEBUG_NAVIGATION
  if( fVerbose > 1 )
  {
    G4cout << "   ---- Exiting ComputeSafety  " << G4endl;
    if( fVerbose > 2 )  { PrintState(); }
    G4cout << "    Returned value of Safety = " << newSafety << G4endl;
  }
  G4cout.precision(oldcoutPrec);
#endif

  return newSafety;
}


// ********************************************************************
// RecheckDistanceToCurrentBoundary
//
// Trial method for checking potential displacement for MS
// Check position aDisplacedGlobalpoint, to see whether it is in the 
// current volume (mother).
// If in mother, check distance to boundary along aNewDirection.
// If in entering daughter, check distance back to boundary. 
// NOTE:
// Can be called only after ComputeStep() is called - before ReLocation
// Deals only with current volume (and potentially entered)
// Caveats
// First VERSION: Does not consider daughter volumes if it remained in mother
//       neither for checking whether it has exited current (mother) volume,
//       nor for determining distance to exit this (mother) volume.
// ********************************************************************
//
G4bool G4Navigator::RecheckDistanceToCurrentBoundary(
                     const G4ThreeVector &aDisplacedGlobalPoint,
                     const G4ThreeVector &aNewDirection,
                     const G4double ProposedMove,
                     G4double *prDistance,
                     G4double *prNewSafety) const
{
  G4ThreeVector localPosition  = ComputeLocalPoint(aDisplacedGlobalPoint);
  G4ThreeVector localDirection = ComputeLocalAxis(aNewDirection);
  // G4double Step = kInfinity;

  G4bool validExitNormal;
  G4ThreeVector exitNormal;
  // Check against mother solid
  G4VPhysicalVolume  *motherPhysical = fHistory.GetTopVolume();
  G4LogicalVolume *motherLogical = motherPhysical->GetLogicalVolume();
  
#ifdef CHECK_ORDER_OF_METHODS
  if( ! fLastTriedStepComputation )
  {
     G4Exception("G4Navigator::RecheckDistanceToCurrentBoundary()",
                 "GeomNav0001", FatalException, 
     "Method must be called after ComputeStep(), before call to LocateMethod.");
  }
#endif

  EInside locatedDaug; // = kUndefined;
  G4double daughterStep= DBL_MAX;
  G4double daughterSafety= DBL_MAX;

  if( fEnteredDaughter )
  {
     if( motherLogical->CharacteriseDaughters() == kReplica )  { return false; }

     // Track arrived at boundary of a daughter volume at 
     //   the last call of ComputeStep().
     // In case the proposed displaced point is inside this daughter,
     //   it must backtrack at least to the entry point.
     // NOTE: No check is made against other daughter volumes.  It is 
     //   assumed that the proposed displacement is small enough that 
     //   this is not needed.

     // Must check boundary of current daughter
     //
     G4VPhysicalVolume *candPhysical= fBlockedPhysicalVolume; 
     G4LogicalVolume *candLogical= candPhysical->GetLogicalVolume();
     G4VSolid        *candSolid=   candLogical->GetSolid();

     G4AffineTransform nextLevelTrf(candPhysical->GetRotation(),
                                    candPhysical->GetTranslation());

     G4ThreeVector dgPosition=  nextLevelTrf.TransformPoint(localPosition); 
     G4ThreeVector dgDirection= nextLevelTrf.TransformAxis(localDirection);
     locatedDaug = candSolid->Inside(dgPosition);

     if( locatedDaug == kInside )
     {
        // Reverse direction - and find first exit. ( Is it valid?)
        // Must backtrack
        G4double distanceBackOut = 
          candSolid->DistanceToOut(dgPosition,
                                   - dgDirection,  // Reverse direction
                                   true, &validExitNormal, &exitNormal);
        daughterStep= - distanceBackOut;
          // No check is made whether the particle could have arrived at 
          // at this location without encountering another volume or 
          // a different psurface of the current volume
        if( prNewSafety )
        {
           daughterSafety= candSolid->DistanceToOut(dgPosition);
        }
     }
     else
     {
        if( locatedDaug == kOutside )
        {
            // See whether it still intersects it
            //
            daughterStep=  candSolid->DistanceToIn(dgPosition,
                                                   dgDirection);
           if( prNewSafety )
           {
              daughterSafety= candSolid->DistanceToIn(dgPosition);
           }
        }
        else
        {
           // The point remains on the surface of candidate solid
           //
           daughterStep= 0.0;
           daughterSafety= 0.0; 
        }
     }

     //  If trial point is in daughter (or on its surface) we have the
     //  answer, the rest is not relevant
     //
     if( locatedDaug != kOutside )
     {
        *prDistance= daughterStep;
        if( prNewSafety )  { *prNewSafety= daughterSafety; }
        return true;
     }
     // If ever extended, so that some type of mother cut daughter, 
     // this would change
  }

  G4VSolid *motherSolid= motherLogical->GetSolid();

  G4double motherStep= DBL_MAX, motherSafety= DBL_MAX;
  
  // Check distance to boundary of mother
  //
  if ( fHistory.GetTopVolumeType()!=kReplica )
  {
     EInside locatedMoth = motherSolid->Inside(localPosition);

     if( locatedMoth == kInside )
     {
        motherSafety= motherSolid->DistanceToOut(localPosition);
        if( ProposedMove >= motherSafety )
        {
           motherStep= motherSolid->DistanceToOut(localPosition,
                                                  localDirection,
                             true, &validExitNormal, &exitNormal);
        }
        else
        {
           motherStep= ProposedMove;
        }
     }
     else if( locatedMoth == kOutside)
     {
        motherSafety= motherSolid->DistanceToIn(localPosition);
        if( ProposedMove >= motherSafety )
        {
            motherStep= - motherSolid->DistanceToIn(localPosition,
                                                   -localDirection);
        }
     }
     else
     {
        motherSafety= 0.0; 
        *prDistance= 0.0;  //  On surface - no move // motherStep;
        if( prNewSafety )  { *prNewSafety= motherSafety; }
        return false;
     }
  }
  else
  {
     return false;
  }
 
  *prDistance=  std::min( motherStep, daughterStep ); 
  if( prNewSafety )
  {
     *prNewSafety= std::min( motherSafety, daughterSafety );
  }
  return true;
}


// ********************************************************************
// CreateTouchableHistoryHandle
// ********************************************************************
//
G4TouchableHistoryHandle G4Navigator::CreateTouchableHistoryHandle() const
{
  return G4TouchableHistoryHandle( CreateTouchableHistory() );
}

// ********************************************************************
// PrintState
// ********************************************************************
//
void  G4Navigator::PrintState() const
{
  G4int oldcoutPrec = G4cout.precision(4);
  if( fVerbose >= 4 )
  {
    G4cout << "The current state of G4Navigator is: " << G4endl;
    G4cout << "  ValidExitNormal= " << fValidExitNormal // << G4endl
           << "  ExitNormal     = " << fExitNormal      // << G4endl
           << "  Exiting        = " << fExiting         // << G4endl
           << "  Entering       = " << fEntering        // << G4endl
           << "  BlockedPhysicalVolume= " ;
    if (fBlockedPhysicalVolume==0)
    {
      G4cout << "None";
    }
    else
    {
      G4cout << fBlockedPhysicalVolume->GetName();
    }
    G4cout << G4endl
           << "  BlockedReplicaNo     = " <<  fBlockedReplicaNo   //  << G4endl
           << "  LastStepWasZero      = " <<   fLastStepWasZero   //  << G4endl
           << G4endl;   
  }
  if( ( 1 < fVerbose) && (fVerbose < 4) )
  {
    G4cout << G4endl; // Make sure to line up
    G4cout << std::setw(30) << " ExitNormal "  << " "
           << std::setw( 5) << " Valid "       << " "     
           << std::setw( 9) << " Exiting "     << " "      
           << std::setw( 9) << " Entering"     << " " 
           << std::setw(15) << " Blocked:Volume "  << " "   
           << std::setw( 9) << " ReplicaNo"        << " "  
           << std::setw( 8) << " LastStepZero  "   << " "   
           << G4endl;   
    G4cout << "( " << std::setw(7) << fExitNormal.x() 
           << ", " << std::setw(7) << fExitNormal.y()
           << ", " << std::setw(7) << fExitNormal.z() << " ) "
           << std::setw( 5)  << fValidExitNormal  << " "   
           << std::setw( 9)  << fExiting          << " "
           << std::setw( 9)  << fEntering         << " ";
    if ( fBlockedPhysicalVolume==0 )
    { G4cout << std::setw(15) << "None"; }
    else
    { G4cout << std::setw(15)<< fBlockedPhysicalVolume->GetName(); }
    G4cout << std::setw( 9)  << fBlockedReplicaNo  << " "
           << std::setw( 8)  << fLastStepWasZero   << " "
           << G4endl;   
  }
  if( fVerbose > 2 ) 
  {
    G4cout.precision(8);
    G4cout << " Current Localpoint = " << fLastLocatedPointLocal << G4endl;
    G4cout << " PreviousSftOrigin  = " << fPreviousSftOrigin << G4endl;
    G4cout << " PreviousSafety     = " << fPreviousSafety << G4endl; 
  }
  G4cout.precision(oldcoutPrec);
}

// ********************************************************************
// ComputeStepLog
// ********************************************************************
//
void G4Navigator::ComputeStepLog(const G4ThreeVector& pGlobalpoint,
                                       G4double moveLenSq) const
{
  //  The following checks only make sense if the move is larger
  //  than the tolerance.

  const G4double fAccuracyForWarning   = kCarTolerance,
                 fAccuracyForException = 1000*kCarTolerance;

  G4ThreeVector OriginalGlobalpoint = fHistory.GetTopTransform().Inverse().
                                      TransformPoint(fLastLocatedPointLocal); 

  G4double shiftOriginSafSq = (fPreviousSftOrigin-pGlobalpoint).mag2();

  // Check that the starting point of this step is 
  // within the isotropic safety sphere of the last point
  // to a accuracy/precision  given by fAccuracyForWarning.
  //   If so give warning.
  //   If it fails by more than fAccuracyForException exit with error.
  //
  if( shiftOriginSafSq >= sqr(fPreviousSafety) )
  {
    G4double shiftOrigin = std::sqrt(shiftOriginSafSq);
    G4double diffShiftSaf = shiftOrigin - fPreviousSafety;

    if( diffShiftSaf > fAccuracyForWarning )
    {
      G4int oldcoutPrec= G4cout.precision(8);
      G4int oldcerrPrec= G4cerr.precision(10);
      std::ostringstream message, suggestion;
      message << "Accuracy error or slightly inaccurate position shift."
              << G4endl
              << "     The Step's starting point has moved " 
              << std::sqrt(moveLenSq)/mm << " mm " << G4endl
              << "     since the last call to a Locate method." << G4endl
              << "     This has resulted in moving " 
              << shiftOrigin/mm << " mm " 
              << " from the last point at which the safety " 
              << "     was calculated " << G4endl
              << "     which is more than the computed safety= " 
              << fPreviousSafety/mm << " mm  at that point." << G4endl
              << "     This difference is " 
              << diffShiftSaf/mm << " mm." << G4endl
              << "     The tolerated accuracy is "
              << fAccuracyForException/mm << " mm.";

      suggestion << " ";
      static G4ThreadLocal G4int warnNow = 0;
      if( ((++warnNow % 100) == 1) )
      {
        message << G4endl
               << "  This problem can be due to either " << G4endl
               << "    - a process that has proposed a displacement"
               << " larger than the current safety , or" << G4endl
               << "    - inaccuracy in the computation of the safety";
        suggestion << "We suggest that you " << G4endl
                   << "   - find i) what particle is being tracked, and "
                   << " ii) through what part of your geometry " << G4endl
                   << "      for example by re-running this event with "
                   << G4endl
                   << "         /tracking/verbose 1 "  << G4endl
                   << "    - check which processes you declare for"
                   << " this particle (and look at non-standard ones)"
                   << G4endl
                   << "   - in case, create a detailed logfile"
                   << " of this event using:" << G4endl
                   << "         /tracking/verbose 6 ";
      }
      G4Exception("G4Navigator::ComputeStep()",
                  "GeomNav1002", JustWarning,
                  message, G4String(suggestion.str()));
      G4cout.precision(oldcoutPrec);
      G4cerr.precision(oldcerrPrec);
    }
#ifdef G4DEBUG_NAVIGATION
    else
    {
      G4cerr << "WARNING - G4Navigator::ComputeStep()" << G4endl
             << "          The Step's starting point has moved "
             << std::sqrt(moveLenSq) << "," << G4endl
             << "          which has taken it to the limit of"
             << " the current safety. " << G4endl;
    }
#endif
  }
  G4double safetyPlus = fPreviousSafety + fAccuracyForException;
  if ( shiftOriginSafSq > sqr(safetyPlus) )
  {
    std::ostringstream message;
    message << "May lead to a crash or unreliable results." << G4endl
            << "        Position has shifted considerably without"
            << " notifying the navigator !" << G4endl
            << "        Tolerated safety: " << safetyPlus << G4endl
            << "        Computed shift  : " << shiftOriginSafSq;
    G4Exception("G4Navigator::ComputeStep()", "GeomNav1002",
                JustWarning, message);
  }
}

// ********************************************************************
// Operator <<
// ********************************************************************
//
std::ostream& operator << (std::ostream &os,const G4Navigator &n)
{
  //  Old version did only the following:
  // os << "Current History: " << G4endl << n.fHistory;
  //  Old behaviour is recovered for fVerbose = 0
  
  // Adapted from G4Navigator::PrintState() const

  G4int oldcoutPrec = os.precision(4);
  if( n.fVerbose >= 4 )
  {
    os << "The current state of G4Navigator is: " << G4endl;
    os << "  ValidExitNormal= " << n.fValidExitNormal << G4endl
    << "  ExitNormal     = " << n.fExitNormal      << G4endl
    << "  Exiting        = " << n.fExiting         << G4endl
    << "  Entering       = " << n.fEntering        << G4endl
    << "  BlockedPhysicalVolume= " ;
    if (n.fBlockedPhysicalVolume==0)
      os << "None";
    else
      os << n.fBlockedPhysicalVolume->GetName();
    os << G4endl
    << "  BlockedReplicaNo     = " <<  n.fBlockedReplicaNo       << G4endl
    << "  LastStepWasZero      = " <<   n.fLastStepWasZero       << G4endl
    << G4endl;
  }
  if( ( 1 < n.fVerbose) && (n.fVerbose < 4) )
  {
    os << G4endl; // Make sure to line up
    os << std::setw(30) << " ExitNormal "  << " "
    << std::setw( 5) << " Valid "       << " "
    << std::setw( 9) << " Exiting "     << " "
    << std::setw( 9) << " Entering"     << " "
    << std::setw(15) << " Blocked:Volume "  << " "
    << std::setw( 9) << " ReplicaNo"        << " "
    << std::setw( 8) << " LastStepZero  "   << " "
    << G4endl;
    os << "( " << std::setw(7) << n.fExitNormal.x()
    << ", " << std::setw(7) << n.fExitNormal.y()
    << ", " << std::setw(7) << n.fExitNormal.z() << " ) "
    << std::setw( 5)  << n.fValidExitNormal  << " "
    << std::setw( 9)  << n.fExiting          << " "
    << std::setw( 9)  << n.fEntering         << " ";
    if ( n.fBlockedPhysicalVolume==0 )
      { os << std::setw(15) << "None"; }
    else
      { os << std::setw(15)<< n.fBlockedPhysicalVolume->GetName(); }
    os << std::setw( 9)  << n.fBlockedReplicaNo  << " "
    << std::setw( 8)  << n.fLastStepWasZero   << " "
    << G4endl;
  }
  if( n.fVerbose > 2 )
  {
    os.precision(8);
    os << " Current Localpoint = " << n.fLastLocatedPointLocal << G4endl;
    os << " PreviousSftOrigin  = " << n.fPreviousSftOrigin << G4endl;
    os << " PreviousSafety     = " << n.fPreviousSafety << G4endl;
  }
  if( n.fVerbose > 3 || n.fVerbose == 0 )
  {
    os << "Current History: " << G4endl << n.fHistory;
  }
    
  os.precision(oldcoutPrec);
  return os;
}
