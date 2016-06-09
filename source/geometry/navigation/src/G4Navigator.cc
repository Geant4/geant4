//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4Navigator.cc,v 1.13 2004/06/18 12:47:05 gcosmo Exp $
// GEANT4 tag $ Name:  $
// 
// class G4Navigator Implementation
//
// Original author: Paul Kent, July 95/96
//
// --------------------------------------------------------------------

#include "G4Navigator.hh"
#include "G4ios.hh"
#include <iomanip>

#include "G4VPhysicalVolume.hh"

// ********************************************************************
// Constructor
// ********************************************************************
//
G4Navigator::G4Navigator()
  : fWasLimitedByGeometry(false), fTopPhysical(0),
    fCheck(false), fPushed(false), fVerbose(0)
{
  ResetStackAndState();

  fActionThreshold_NoZeroSteps  = 10; 
  fAbandonThreshold_NoZeroSteps = 25; 
}

// ********************************************************************
// Destructor
// ********************************************************************
//
G4Navigator::~G4Navigator()
{;}

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
  
  if( considerDirection && pGlobalDirection != 0 )
  {
    globalDirection=*pGlobalDirection;
  }

#ifdef G4DEBUG_NAVIGATION
  G4cerr << "Upon entering LocateGlobalPointAndSetup():" << G4endl;
  G4cerr << "    History = " << G4endl << fHistory << G4endl << G4endl;
#endif

#ifdef G4VERBOSE
  G4int oldcoutPrec = G4cout.precision(8);
  if( fVerbose > 2 )
  {
    G4cout << "*** G4Navigator::LocateGlobalPointAndSetup: ***" << G4endl; 
    G4cout << "    Called with arguments: " << G4endl
           << "    Globalpoint = " << globalPoint << G4endl
           << "    RelativeSearch = " << relativeSearch  << G4endl;
    if( fVerbose == 4 )
    {
      G4cout << "    ----- Upon entering:" << G4endl;
      PrintState();
    }
  }
#endif

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
        if ( fHistory.GetDepth() )
        {
          fBlockedPhysicalVolume = fHistory.GetTopVolume();
          fBlockedReplicaNo = fHistory.GetTopReplicaNo();
          fHistory.BackLevel();
        }
        else
        {
          fLastLocatedPointLocal = localPoint;
          return 0;           // Have exited world volume
        }
        // A fix for the case where a volume is "entered" at an edge
        // and a coincident surface exists outside it.
        //  - This stops it from exiting further volumes and cycling
        //  - However ReplicaNavigator treats this case itself
        //
        if ( fLocatedOnEdge && (VolumeType(fBlockedPhysicalVolume)!=kReplica ))
        { 
          fExiting= false;
        }
      }
      else
        if ( fEntering )
        {
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
              G4VSolid *pSolid;
              G4VPVParameterisation *pParam;
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
              pLogical->SetMaterial(pParam->ComputeMaterial(fBlockedReplicaNo, 
                                                      fBlockedPhysicalVolume));
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
  while (notKnownContained)
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
      if ( fHistory.GetDepth() )
      {
        fBlockedPhysicalVolume = fHistory.GetTopVolume();
        fBlockedReplicaNo = fHistory.GetTopReplicaNo();
        fHistory.BackLevel();
        fExiting = false;
      }
      else
      {
        fLastLocatedPointLocal = localPoint;
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
          if ( fHistory.GetTopVolumeType()!=kReplica )
          {
            G4ThreeVector normal = targetSolid->SurfaceNormal(localPoint);
            directionExiting = normal.dot(localDirection) > 0.0;
            isExiting = isExiting || directionExiting;
          }
        }
        if( isExiting )
        {
          if ( fHistory.GetDepth() )
          {
            fBlockedPhysicalVolume = fHistory.GetTopVolume();
            fBlockedReplicaNo = fHistory.GetTopReplicaNo();
            fHistory.BackLevel();
            //
            // Still on surface but exited volume not necessarily convex
            //
            fValidExitNormal = false;
          } 
          else
          {
            fLastLocatedPointLocal = localPoint;
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
        noResult = fparamNav.LevelLocate(fHistory,
                                         fBlockedPhysicalVolume,
                                         fBlockedReplicaNo,
                                         globalPoint,
                                         pGlobalDirection,
                                         considerDirection,
                                         localPoint);
        break;
    }

    // LevelLocate returns true if it finds a daughter volume 
    // in which globalPoint is inside (or on the surface).

    if ( noResult )
    {
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
#ifdef G4VERBOSE
      if( fVerbose > 1 )
      { 
         G4VPhysicalVolume* enteredPhysical = fHistory.GetTopVolume();
         G4cout << "*** G4Navigator::LocateGlobalPointAndSetup: ***" << G4endl; 
         G4cout << "    Entering volume: " << enteredPhysical->GetName()
                << G4endl;
      }
#endif
    }
  } while (noResult);

  fLastLocatedPointLocal = localPoint;

#ifdef G4VERBOSE
  if( fVerbose == 4 )
  {
    G4cout.precision(6);
    G4String curPhysVol_Name("None");
    if (targetPhysical!=0)
      curPhysVol_Name = targetPhysical->GetName();
    G4cout << "    Return value = new volume = " << curPhysVol_Name << G4endl;
    G4cout << "    ----- Upon exiting:" << G4endl;
    PrintState();
#ifdef G4DEBUG_NAVIGATION
    G4cerr << "Upon exiting LocateGlobalPointAndSetup():" << G4endl;
    G4cerr << "    History = " << G4endl << fHistory << G4endl << G4endl;
#endif
  }
  G4cout.precision(oldcoutPrec);
#endif

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
   fLastLocatedPointLocal = ComputeLocalPoint(pGlobalpoint);

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
         //  else { fnormalNav. nothing !? }
         break;

       case kParameterised:
         // Resets state & returns voxel node
         //
         fparamNav.ParamVoxelLocate( pVoxelHeader, fLastLocatedPointLocal );
         break;

       case kReplica:
         G4Exception("G4Navigator::LocateGlobalPointWithinVolume()",
                     "NotApplicable", FatalException,
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
// LocateGlobalPointAndUpdateTouchableHandle
// ********************************************************************
//
void G4Navigator::LocateGlobalPointAndUpdateTouchableHandle(
                               const G4ThreeVector&       position,
                               const G4ThreeVector&       direction,
                                     G4TouchableHandle&   oldTouchableToUpdate,
                               const G4bool               RelativeSearch )
{
  G4VPhysicalVolume* pPhysVol;
  pPhysVol = LocateGlobalPointAndSetup( position,&direction,RelativeSearch );
  if( fEnteredDaughter || fExitedMother )
  {
     oldTouchableToUpdate = CreateTouchableHistory();
     if( pPhysVol == 0 )
     {
       // We want to ensure that the touchable is correct in this case.
       //  The method below should do this and recalculate a lot more ....
       //
       oldTouchableToUpdate->UpdateYourself( pPhysVol, &fHistory );
     }
  }
  return;
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
// fLastStepWasZero  - True if last Step size was zero.
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
// fLastStepWasZero  - True if this Step size was zero.
// ********************************************************************
//
G4double G4Navigator::ComputeStep( const G4ThreeVector &pGlobalpoint,
                                   const G4ThreeVector &pDirection,
                                   const G4double pCurrentProposedStepLength,
                                         G4double &pNewSafety)
{
  G4ThreeVector localDirection = ComputeLocalAxis(pDirection);
  G4double Step = DBL_MAX;
  G4VPhysicalVolume  *motherPhysical = fHistory.GetTopVolume();
  G4LogicalVolume *motherLogical = motherPhysical->GetLogicalVolume();

  static G4int sNavCScalls=0;
  sNavCScalls++;

#ifdef G4VERBOSE
  G4int oldcoutPrec= G4cout.precision(8);
  G4int oldcerrPrec= G4cerr.precision(10);
  if( fVerbose > 0 )
  {
    G4cout << "*** G4Navigator::ComputeStep: ***" << G4endl; 
    G4cout << "    Volume = " << motherPhysical->GetName() 
           << " - Proposed step length = " << pCurrentProposedStepLength
           << G4endl; 
    if( fVerbose == 4 ) 
    {
      G4cout << "    Called with the arguments: " << G4endl
             << "    Globalpoint = " << std::setw(25) << pGlobalpoint
             << G4endl
             << "    Direction   = " << std::setw(25) << pDirection
             << G4endl;
      G4cout << "    ----- Upon entering :" << G4endl;
      PrintState();
    }
  }

  static G4double fAccuracyForWarning   = kCarTolerance,
                  fAccuracyForException = 1000*kCarTolerance;
#endif

  G4ThreeVector newLocalPoint = ComputeLocalPoint(pGlobalpoint);
  if( newLocalPoint != fLastLocatedPointLocal )
  {
    // Check whether the relocation is within safety
    //
    G4ThreeVector oldLocalPoint = fLastLocatedPointLocal;
    G4double moveLenSq = (newLocalPoint-oldLocalPoint).mag2();

    if ( moveLenSq >= kCarTolerance*kCarTolerance )
    {
#ifdef G4VERBOSE
      //  The following checks only make sense if the move is larger
      //  than the tolerance.
      // 
      G4ThreeVector OriginalGlobalpoint =
                    fHistory.GetTopTransform().Inverse().
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
        G4double shiftOrigin = sqrt(shiftOriginSafSq);
        G4double diffShiftSaf = shiftOrigin - fPreviousSafety;

        if( diffShiftSaf > fAccuracyForWarning )
        {
          G4Exception("G4Navigator::ComputeStep()",
                      "UnexpectedPositionShift", JustWarning,
                      "Accuracy ERROR or slightly inaccurate position shift.");
          G4cerr << "     The Step's starting point has moved " 
                 << sqrt(moveLenSq)/mm << " mm " << G4endl
                 << "     since the last call to a Locate method." << G4endl;
          G4cerr << "     This has resulted in moving " 
                 << shiftOrigin/mm << " mm " 
                 << " from the last point at which the safety " 
                 << "     was calculated " << G4endl;
          G4cerr << "     which is more than the computed safety= " 
                 << fPreviousSafety/mm << " mm  at that point." << G4endl;
          G4cerr << "     This difference is " 
                 << diffShiftSaf/mm << " mm." << G4endl
                 << "     The tolerated accuracy is "
                 << fAccuracyForException/mm << " mm." << G4endl;

          static G4int warnNow = 0;
          if( ((++warnNow % 100) == 1) )
          {
            G4cerr << "  This problem can be due to either " << G4endl;
            G4cerr << "    - a process that has proposed a displacement"
                   << " larger than the current safety , or" << G4endl;
            G4cerr << "    - inaccuracy in the computation of the safety"
                   << G4endl;
            G4cerr << "  We suggest that you " << G4endl
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
                   << "         /tracking/verbose 6 "
                   << G4endl;
          }
        }
#ifdef G4DEBUG_NAVIGATION
        else
        {
          G4cerr << "WARNING - G4Navigator::ComputeStep()" << G4endl
                 << "          The Step's starting point has moved "
                 << sqrt(moveLenSq) << "," << G4endl
                 << "          which has taken it to the limit of"
                 << " the current safety. " << G4endl;
        }
#endif
      }
      G4double safetyPlus = fPreviousSafety + fAccuracyForException;
      if ( shiftOriginSafSq > sqr(safetyPlus) )
      {
        G4cerr << "ERROR - G4Navigator::ComputeStep()" << G4endl
               << "        Position has shifted considerably without"
               << " notifying the navigator !" << G4endl
               << "        Tolerated safety: " << safetyPlus << G4endl
               << "        Computed shift  : " << shiftOriginSafSq << G4endl;
        G4Exception("G4Navigator::ComputeStep()",
                    "SignificantPositionShift", JustWarning,
                    "May lead to a crash or unreliable results.");
      }
#endif  // end G4VERBOSE

      // Relocate the point within the same volume
      //
      LocateGlobalPointWithinVolume( pGlobalpoint );
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
        break;
      case kParameterised:
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
        break;
      case kReplica:
        G4Exception("G4Navigator::ComputeStep()", "NotApplicable",
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
    Step = freplicaNav.ComputeStep(pGlobalpoint,
                                   pDirection,
                                   fLastLocatedPointLocal,
                                   localDirection,
                                   pCurrentProposedStepLength,
                                   pNewSafety,
                                   fHistory,
                                   fValidExitNormal,
                                   fExitNormal,
                                   exitingReplica,
                                   fEntering,
                                   &fBlockedPhysicalVolume,
                                   fBlockedReplicaNo);
    fExiting= exitingReplica;                          // still ok to set it ??
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
  fLastStepWasZero = (Step==0.0);
  if (fPushed)  fPushed = fLastStepWasZero;

  // Handle large number of consecutive zero steps
  //
  if ( fLastStepWasZero )
  {
    fNumberZeroSteps++;
#ifdef G4DEBUG_NAVIGATION
    if( fNumberZeroSteps > 1 )
    {
       G4cout << "G4Nav - CompStep: another zero step, # " << fNumberZeroSteps
              << " at " << pGlobalpoint
              << " in volume " << motherPhysical->GetName()
              << " nav-comp-step calls # " << sNavCScalls
              << G4endl;
    }
#endif
    if( (fNumberZeroSteps > fActionThreshold_NoZeroSteps-1) && (!fPushed) )
    {
       // Act to recover this stuck track. Pushing it once along direction
       //
       Step += 0.9*kCarTolerance;
       fPushed = true;
#ifdef G4VERBOSE
       G4cerr << "WARNING - G4Navigator::ComputeStep()" << G4endl
              << "          Track stuck, not moving for " 
              << fNumberZeroSteps << " steps" << G4endl
              << "          in volume -" << motherPhysical->GetName()
              << "- at point " << pGlobalpoint << G4endl
              << "          direction: " << pDirection << "." << G4endl
              << "          Potential geometry or navigation problem !"
              << G4endl
              << "          Trying pushing it of " << Step << " mm ..."
              << G4endl;
#endif
     }
     if( fNumberZeroSteps > fAbandonThreshold_NoZeroSteps-1 )
     {
        // Must kill this stuck track
        //
        G4cerr << "ERROR - G4Navigator::ComputeStep()" << G4endl
               << "        Track stuck, not moving for " 
               << fNumberZeroSteps << " steps" << G4endl
               << "        in volume -" << motherPhysical->GetName()
               << "- at point " << pGlobalpoint << G4endl
               << "        direction: " << pDirection << "." << G4endl;
        G4Exception("G4Navigator::ComputeStep()",
                    "StuckTrack", EventMustBeAborted, 
                    "Stuck Track: potential geometry or navigation problem.");
    }
  }
  else
  {
    if (!fPushed)  fNumberZeroSteps = 0;
  }

  fEnteredDaughter = fEntering;   // I expect to enter a volume in this Step
  fExitedMother = fExiting;

  if( fExiting )
  {
#ifdef G4DEBUG_NAVIGATION
    G4cout << " At G4Nav CompStep End - if(exiting) - fExiting= " << fExiting 
           << " fValidExitNormal = " << fValidExitNormal  << G4endl;
    G4cout << " fExitNormal= " << fExitNormal << G4endl;
#endif

    if(fValidExitNormal)
    {
      // Convention: fExitNormal is in the 'grand-mother' coordinate system
      //
      fGrandMotherExitNormal= fExitNormal;

      // If no relocation were made, we would need to rotate it back to 
      // this (the mother) coordinate system
      // const G4RotationMatrix* motherRotation= motherPhysical->GetRotation();
      // G4ThreeVector trueMotherExitNormal = fGrandMotherExitNormal;
      // Un-rotate gran->mother
      //
      // trueMotherExitNormal *= (*motherRotation);

      // However, relocation will put us either in
      //   - the grand-mother (OK)
      //   - in the grand-grand mother (to BE checked if this is dealt with)
    }
    else
    {  
      // We must calculate the normal anyway
      // (in order to have it if requested)
      //
      G4ThreeVector finalGlobalPoint, finalLocalPoint, localExitNormal;
      finalLocalPoint = fLastLocatedPointLocal + localDirection*Step;
      localExitNormal  = motherLogical->GetSolid()->
                         SurfaceNormal(finalLocalPoint);
      const G4RotationMatrix* mRot = motherPhysical->GetRotation();
      if( mRot )
      { 
         G4ThreeVector grandMotherExitNormal = localExitNormal;
         grandMotherExitNormal *= (*mRot);

         // Now fGrandMotherExitNormal is in the 'grand-mother'
         // coordinate system
         //
         fGrandMotherExitNormal = grandMotherExitNormal;
      }
    }
#ifdef G4DEBUG_NAVIGATION
    G4cout << " fGrandMotherExitNormal= " << fGrandMotherExitNormal << G4endl;
#endif
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
    if( fVerbose == 4 ) 
    {
      G4cout << "    ----- Upon exiting :" << G4endl;
      PrintState();
    }
    G4cout <<"    Returned step = " << Step << G4endl;
    if( Step == kInfinity )
      G4cout << "    Original proposed step = "
             << pCurrentProposedStepLength << G4endl;
    G4cout << "    Safety = " << pNewSafety << G4endl;
  }
  G4cout.precision(oldcoutPrec);
  G4cerr.precision(oldcerrPrec);
#endif

  return Step;
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
  fExitNormal            = G4ThreeVector(0,0,0);

  fPreviousSftOrigin     = G4ThreeVector(0,0,0);
  fPreviousSafety        = 0.0; 

  fNumberZeroSteps       = 0;
    
  fBlockedPhysicalVolume = 0;
  fBlockedReplicaNo      = -1;
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
  G4VPhysicalVolume *mother, *current;
  G4VSolid *pSolid;
  G4VPVParameterisation *pParam;

  mother = fHistory.GetVolume(0);
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

        // Set up the correct solid and material in Logical Volume
        //
        G4LogicalVolume *pLogical = current->GetLogicalVolume();
        pLogical->SetSolid( pSolid );
        pLogical->SetMaterial( pParam->ComputeMaterial(replicaNo, current));
        break;
    }
    mother = current;
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
  G4ThreeVector ExitNormal(0.,0.,0.);

  if ( EnteredDaughterVolume() )
  {
    ExitNormal= -(fHistory.GetTopVolume()->GetLogicalVolume()->
                  GetSolid()->SurfaceNormal(fLastLocatedPointLocal));
    *valid = true;
  }
  else
    if( fExitedMother )
    {
      ExitNormal = fGrandMotherExitNormal;
      *valid = true;
    }
    else
    {
      // We are not at a boundary.
      // ExitNormal remains (0,0,0)
      //
      *valid = false;
    }
  return ExitNormal;
}

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
                                     const G4double pMaxLength)
{
  G4double newSafety = 0.0;

#ifdef G4VERBOSE
  G4int oldcoutPrec = G4cout.precision(8);
  if( fVerbose > 0 )
  {
    G4VPhysicalVolume  *motherPhysical = fHistory.GetTopVolume();
    G4cout << "*** G4Navigator::ComputeSafety: ***" << G4endl; 
    G4cout << "    Volume = " << motherPhysical->GetName() 
           << " - Maximum length = " << pMaxLength << G4endl; 
    if( fVerbose == 4 ) 
    {
      G4cout << "    Called with the arguments: " << G4endl
             << "    Globalpoint = " << pGlobalpoint << G4endl;
      G4cout << "    ----- Upon entering :" << G4endl;
      PrintState();
    }
  }
#endif

  if( !(fEnteredDaughter || fExitedMother) )
  {
    // Pseudo-relocate to this point (updates voxel information only)
    //
    LocateGlobalPointWithinVolume( pGlobalpoint );

    G4VPhysicalVolume *motherPhysical = fHistory.GetTopVolume();
    G4LogicalVolume *motherLogical = motherPhysical->GetLogicalVolume();
    G4ThreeVector localPoint = ComputeLocalPoint(pGlobalpoint);

    if ( fHistory.GetTopVolumeType()!=kReplica )
    {
      switch(CharacteriseDaughters(motherLogical))
      {
        case kNormal:
          if ( motherLogical->GetVoxelHeader() )
          {
            newSafety=fvoxelNav.ComputeSafety(localPoint,fHistory,pMaxLength);
          }
          else
          {
            newSafety=fnormalNav.ComputeSafety(localPoint,fHistory,pMaxLength);
          }
          break;
        case kParameterised:
          newSafety = fparamNav.ComputeSafety(localPoint,fHistory,pMaxLength);
          break;
        case kReplica:
          G4Exception("G4Navigator::ComputeSafety()", "NotApplicable",
                      FatalException, "Not applicable for replicated volumes.");
          break;
      }
    }
    else
    {
      newSafety = freplicaNav.ComputeSafety(pGlobalpoint, localPoint,
                                            fHistory, pMaxLength);
    }
  }

  // Remember last safety origin & value
  //
  fPreviousSftOrigin = pGlobalpoint;
  fPreviousSafety = newSafety; 

#ifdef G4VERBOSE
  if( fVerbose > 1 ) 
  {
    G4cout << "    ----- Upon exiting :" << G4endl;
    PrintState();
    G4cout << "    Returned value of Safety = " << newSafety << G4endl;
  }
  G4cout.precision(oldcoutPrec);
#endif

  return newSafety;
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
void  G4Navigator::PrintState()
{
  G4int oldcoutPrec = G4cout.precision(4);
  if( fVerbose == 4 )
  {
    G4cout << "The current state of G4Navigator is: " << G4endl;
    G4cout << "  ValidExitNormal= " << fValidExitNormal << G4endl
           << "  ExitNormal     = " << fExitNormal      << G4endl
           << "  Exiting        = " << fExiting         << G4endl
           << "  Entering       = " << fEntering        << G4endl
           << "  BlockedPhysicalVolume= " ;
    if (fBlockedPhysicalVolume==0)
      G4cout << "None";
    else
      G4cout << fBlockedPhysicalVolume->GetName();
    G4cout << G4endl
           << "  BlockedReplicaNo     = " <<  fBlockedReplicaNo       << G4endl
           << "  LastStepWasZero      = " <<   fLastStepWasZero       << G4endl
           << G4endl;   
  }
  if( ( 1 < fVerbose) && (fVerbose < 4) )
  {
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
      G4cout << std::setw(15) << "None";
    else
      G4cout << std::setw(15)<< fBlockedPhysicalVolume->GetName();
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
// Operator <<
// ********************************************************************
//
std::ostream& operator << (std::ostream &os,const G4Navigator &n)
{
  os << "Current History: " << G4endl << n.fHistory;
  return os;
}
