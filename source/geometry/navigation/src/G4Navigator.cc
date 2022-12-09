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
// class G4Navigator Implementation
//
// Original author: Paul Kent, July 95/96
// Responsible 1996-present: John Apostolakis, Gabriele Cosmo
// Additional revisions by: Pedro Arce, Vladimir Grichine
// --------------------------------------------------------------------

#include <iomanip>

#include "G4Navigator.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4GeometryTolerance.hh"
#include "G4VPhysicalVolume.hh"

#include "G4VoxelSafety.hh"

// Constant determining how precise normals should be (how close to unit
// vectors). If exceeded, warnings will be issued.
// Can be CLHEP::perMillion (its old default) for geometry checking.
//
static const G4double kToleranceNormalCheck = CLHEP::perThousand;

// ********************************************************************
// Constructor
// ********************************************************************
//
G4Navigator::G4Navigator()
{
  ResetStackAndState();
    // Initialises also all 
    // - exit / entry flags
    // - flags & variables for exit normals
    // - zero step counters
    // - blocked volume 

  if( fVerbose > 2 )
  {
    G4cout << " G4Navigator parameters: Action Threshold (No Zero Steps) = "
           << fActionThreshold_NoZeroSteps
           << "  Abandon Threshold (No Zero Steps) = "
           << fAbandonThreshold_NoZeroSteps << G4endl;
  }
  kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  fMinStep = 0.05*kCarTolerance;
  fSqTol = sqr(kCarTolerance);

  fregularNav.SetNormalNavigation( &fnormalNav );

  fStepEndPoint = G4ThreeVector( kInfinity, kInfinity, kInfinity ); 
  fLastStepEndPointLocal = G4ThreeVector( kInfinity, kInfinity, kInfinity ); 

  fpVoxelSafety = new G4VoxelSafety();
#ifdef ALTERNATIVE_VOXEL_NAV
  fpvoxelNav    = new G4VoxelNavigation();
#endif  
}

// ********************************************************************
// Destructor
// ********************************************************************
//
G4Navigator::~G4Navigator()
{
  delete fpVoxelSafety;
  delete fpExternalNav;
#ifdef ALTERNATIVE_VOXEL_NAV  
  delete fpvoxelNav;
#endif
}

// ********************************************************************
// ResetHierarchyAndLocate
// ********************************************************************
//
G4VPhysicalVolume*
G4Navigator::ResetHierarchyAndLocate(const G4ThreeVector& p,
                                     const G4ThreeVector& direction,
                                     const G4TouchableHistory& h)
{
  ResetState();
  fHistory = *h.GetHistory();
  SetupHierarchy();
  fLastTriedStepComputation = false;  // Redundant, but best
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
  G4bool notKnownContained = true, noResult;
  G4VPhysicalVolume *targetPhysical;
  G4LogicalVolume *targetLogical;
  G4VSolid *targetSolid = 0;
  G4ThreeVector localPoint, globalDirection;
  EInside insideCode;

  G4bool considerDirection = pGlobalDirection && ((!ignoreDirection) || fLocatedOnEdge);

  fLastTriedStepComputation = false;   
  fChangedGrandMotherRefFrame = false;  // For local exit normal
   
  if( considerDirection )
  {
    globalDirection=*pGlobalDirection;
  }

#ifdef G4VERBOSE
  if( fVerbose > 2 )
  {
    G4long oldcoutPrec = G4cout.precision(8);
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

  G4int noLevelsExited = 0;

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
        ++noLevelsExited;  // count this first level entered too

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
          
          return nullptr;           // Have exited world volume
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
          fExiting = false;
          // Consider effect on Exit Normal !?
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
            case kExternal:
              G4Exception("G4Navigator::LocateGlobalPointAndSetup()",
                          "GeomNav0001", FatalException,
                          "Extra levels not applicable for external volumes.");
              break;
          }
          fEntering = false;
          fBlockedPhysicalVolume = nullptr;
          localPoint = fHistory.GetTopTransform().TransformPoint(globalPoint);
          notKnownContained = false;
        }
    }
    else
    {
      fBlockedPhysicalVolume = nullptr;
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
    EVolume topVolumeType = fHistory.GetTopVolumeType(); 
    if (topVolumeType!=kReplica && topVolumeType!=kExternal)
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
       if( topVolumeType == kReplica )
       {
          insideCode = freplicaNav.BackLocate(fHistory, globalPoint, localPoint,
                                              fExiting, notKnownContained);
          // !CARE! if notKnownContained returns false then the point is within
          // the containing placement volume of the replica(s). If insidecode
          // will result in the history being backed up one level, then the
          // local point returned is the point in the system of this new level
       }
       else
       {
          targetSolid = fHistory.GetTopVolume()->GetLogicalVolume()->GetSolid();
          localPoint = fHistory.GetTopTransform().TransformPoint(globalPoint);
          G4ThreeVector localDirection =
             fHistory.GetTopTransform().TransformAxis(globalDirection);
          insideCode = fpExternalNav->Inside(targetSolid, localPoint, localDirection);
       }
    }

    // Point is inside current volume, break out of the loop
    if ( insideCode == kInside )
      break;

    // Point is outside current volume, move up a level in the hierarchy
    if ( insideCode == kOutside )
    {
      ++noLevelsExited;

      // Exiting world volume
      if ( fHistory.GetDepth() == 0 )
      {
        fLocatedOutsideWorld = true;
        fLastLocatedPointLocal = localPoint;
        return nullptr;
      }

      fBlockedPhysicalVolume = fHistory.GetTopVolume();
      fBlockedReplicaNo = fHistory.GetTopReplicaNo();
      fHistory.BackLevel();
      fExiting = false;

      if( noLevelsExited > 1 )
      {
        // The first transformation was done by the sub-navigator
        //
        if(const auto *mRot = fBlockedPhysicalVolume->GetRotation())
        {
          fGrandMotherExitNormal *= (*mRot).inverse();
          fChangedGrandMotherRefFrame = true;
        }
      }
      continue;
    }

    // Point is on the surface of a volume
    G4bool isExiting = fExiting;
    if( (!fExiting) && considerDirection )
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
      if ( fHistory.GetTopVolumeType() != kReplica )
      {
        G4ThreeVector normal = targetSolid->SurfaceNormal(localPoint);
        directionExiting = normal.dot(localDirection) > 0.0;
        isExiting = isExiting || directionExiting;
      }
    }

    // Point is on a surface, but no longer exiting, break out of the loop
    if ( !isExiting )
      break;

    ++noLevelsExited;

    // Point is on the outer surface, leaving world volume
    if ( fHistory.GetDepth() == 0 )
    {
      fLocatedOutsideWorld = true;
      fLastLocatedPointLocal = localPoint;
      return nullptr;
    }

    // Point is still on a surface, but exited a volume not necessarily convex
    fValidExitNormal = false;
    fBlockedPhysicalVolume = fHistory.GetTopVolume();
    fBlockedReplicaNo = fHistory.GetTopReplicaNo();
    fHistory.BackLevel();

    if( noLevelsExited > 1 )
    {
      // The first transformation was done by the sub-navigator
      //
      const G4RotationMatrix* mRot =
        fBlockedPhysicalVolume->GetRotation();
      if( mRot )
      {
        fGrandMotherExitNormal *= (*mRot).inverse();
        fChangedGrandMotherRefFrame = true;
      }
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
          noResult = GetVoxelNavigator().LevelLocate(fHistory,
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
      case kExternal:
        noResult = fpExternalNav->LevelLocate(fHistory,
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
      fBlockedPhysicalVolume = nullptr;
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
    G4long oldcoutPrec = G4cout.precision(8);
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

  fLocatedOutsideWorld = false;

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
   assert( !fWasLimitedByGeometry );   
   // Check: Either step was not limited by a boundary or 
   //          else the full step is no longer being taken
#endif
  
   fLastLocatedPointLocal = ComputeLocalPoint(pGlobalpoint);
   fLastTriedStepComputation = false;
   fChangedGrandMotherRefFrame = false;  //  Frame for Exit Normal

   // For the case of Voxel (or Parameterised) volume the respective 
   // Navigator must be messaged to update its voxel information etc

   // Update the state of the Sub Navigators 
   // - in particular any voxel information they store/cache
   //
   G4VPhysicalVolume*  motherPhysical = fHistory.GetTopVolume();
   G4LogicalVolume*    motherLogical  = motherPhysical->GetLogicalVolume();
   G4SmartVoxelHeader* pVoxelHeader   = motherLogical->GetVoxelHeader();

   switch( CharacteriseDaughters(motherLogical) )
   {
       case kNormal:
         if ( pVoxelHeader )
         {
           GetVoxelNavigator().VoxelLocate( pVoxelHeader, fLastLocatedPointLocal );
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
         // Nothing to do
         break;
       case kExternal:
         fpExternalNav->RelocateWithinVolume( motherPhysical,
                                              fLastLocatedPointLocal );
         break;
   }

   // Reset the state variables 
   //   - which would have been affected
   //     by the 'equivalent' call to LocateGlobalPointAndSetup
   //   - who's values have been invalidated by the 'move'.
   //
   fBlockedPhysicalVolume = nullptr; 
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
  fSaveState.sBlockedReplicaNo = fBlockedReplicaNo;

  fSaveState.sLastStepWasZero = fLastStepWasZero;
  
  fSaveState.sLocatedOutsideWorld = fLocatedOutsideWorld;
  fSaveState.sLastLocatedPointLocal = fLastLocatedPointLocal;
  fSaveState.sEnteredDaughter = fEnteredDaughter;
  fSaveState.sExitedMother = fExitedMother;
  fSaveState.sWasLimitedByGeometry = fWasLimitedByGeometry;

  // Even the safety sphere - if you want to change it do it explicitly!
  //
  fSaveState.sPreviousSftOrigin = fPreviousSftOrigin;
  fSaveState.sPreviousSafety = fPreviousSafety;
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
  fBlockedReplicaNo = fSaveState.sBlockedReplicaNo; 

  fLastStepWasZero = fSaveState.sLastStepWasZero;
  
  fLocatedOutsideWorld = fSaveState.sLocatedOutsideWorld;
  fLastLocatedPointLocal = fSaveState.sLastLocatedPointLocal;
  fEnteredDaughter = fSaveState.sEnteredDaughter;
  fExitedMother = fSaveState.sExitedMother;
  fWasLimitedByGeometry = fSaveState.sWasLimitedByGeometry;
  
  // The 'expected' behaviour is to restore these too (fix 2014.05.26)
  fPreviousSftOrigin = fSaveState.sPreviousSftOrigin;
  fPreviousSafety = fSaveState.sPreviousSafety;
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
G4double G4Navigator::ComputeStep( const G4ThreeVector& pGlobalpoint,
                                   const G4ThreeVector& pDirection,
                                   const G4double pCurrentProposedStepLength,
                                         G4double& pNewSafety)
{
#ifdef G4DEBUG_NAVIGATION
  static G4ThreadLocal G4int sNavCScalls = 0;
  ++sNavCScalls;
#endif

  G4ThreeVector localDirection = ComputeLocalAxis(pDirection);
  G4double Step = kInfinity;
  G4VPhysicalVolume  *motherPhysical = fHistory.GetTopVolume();
  G4LogicalVolume *motherLogical = motherPhysical->GetLogicalVolume();

  // All state relating to exiting normals must be reset
  //
  fExitNormalGlobalFrame = G4ThreeVector( 0., 0., 0.);
    // Reset value - to erase its memory
  fChangedGrandMotherRefFrame = false;
    // Reset - used for local exit normal
  fGrandMotherExitNormal = G4ThreeVector( 0., 0., 0.); 
  fCalculatedExitNormal = false;
    // Reset for new step

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
    }
  }
  if ( fHistory.GetTopVolumeType()!=kReplica )
  {
    switch( CharacteriseDaughters(motherLogical) )
    {
      case kNormal:
        if ( motherLogical->GetVoxelHeader() )
        {
          Step = GetVoxelNavigator().ComputeStep(fLastLocatedPointLocal,
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
      case kExternal:
        Step = fpExternalNav->ComputeStep(fLastLocatedPointLocal,
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
    }
  }
  else
  {
    // In the case of a replica, it must handle the exiting
    // edge/corner problem by itself
    //
    fExiting = fExitedMother;
    Step = freplicaNav.ComputeStep(pGlobalpoint,
                                   pDirection,
                                   fLastLocatedPointLocal,
                                   localDirection,
                                   pCurrentProposedStepLength,
                                   pNewSafety,
                                   fHistory,
                                   fValidExitNormal,
                                   fCalculatedExitNormal,
                                   fExitNormal,
                                   fExiting,
                                   fEntering,
                                   &fBlockedPhysicalVolume,
                                   fBlockedReplicaNo);
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
    ++fNumberZeroSteps;
    
    G4bool act     = fNumberZeroSteps >= fActionThreshold_NoZeroSteps;
    G4bool actAndReport = false;
    G4bool abandon = fNumberZeroSteps >= fAbandonThreshold_NoZeroSteps;
    G4bool inform  = false;
#ifdef G4VERBOSE
    actAndReport = act && (!fPushed) && fWarnPush;
#endif
#ifdef G4DEBUG_NAVIGATION
    inform = fNumberZeroSteps > 1;
#endif
    
    if ( act || inform )
    {
      if( act && !abandon )
      {
        // Act to recover this stuck track. Pushing it along original direction
        //
        Step += 100*kCarTolerance;
        fPushed = true;
      }

      if( actAndReport || abandon || inform )
      {
        std::ostringstream message;
         
        message.precision(16);      
        message << "Stuck Track: potential geometry or navigation problem."
                << G4endl;
        message << "  Track stuck, not moving for " 
                << fNumberZeroSteps << " steps."  << G4endl
                << "  Current  phys volume: '" << motherPhysical->GetName()
                << "'" << G4endl
                << "   - at position : " << pGlobalpoint << G4endl
                << "     in direction: " << pDirection   << G4endl
                << "    (local position: " << newLocalPoint << ")" << G4endl
                << "    (local direction: " << localDirection << ")." << G4endl
                << "  Previous phys volume: '"
                << ( fLastMotherPhys ? fLastMotherPhys->GetName() : "" )
                << "'" << G4endl << G4endl;
        if( actAndReport || abandon )
        {
           message << "  Likely geometry overlap - else navigation problem !"
                   << G4endl;
        }
        if( abandon ) // i.e. fNumberZeroSteps >= fAbandonThreshold_NoZeroSteps
        {
          // Must kill this stuck track
#ifdef G4VERBOSE
          if ( fWarnPush ) { CheckOverlapsIterative(motherPhysical); }
#endif
          message << " Track *abandoned* due to excessive number of Zero steps."
                  << " Event aborted. " << G4endl << G4endl;
          G4Exception("G4Navigator::ComputeStep()", "GeomNav0003",
                      EventMustBeAborted, message);
        }
        else
        {
#ifdef G4VERBOSE
          if ( actAndReport )  // (!fPushed => !wasPushed) && (fWarnPush))
          {
             message << "   *** Trying to get *unstuck* using a push"
                     << " - expanding step to " << Step << " (mm) ..."
                     << "       Potential overlap in geometry !" << G4endl;
             G4Exception("G4Navigator::ComputeStep()", "GeomNav1002",
                         JustWarning, message); 
          }
#endif
#ifdef G4DEBUG_NAVIGATION      
          else
          {
            if( fNumberZeroSteps > 1 )
            {
               message << ", nav-comp-step calls # " << sNavCScalls
                       << ", Step= " << Step << G4endl;
               G4cout << message.str();
            }
          }
#endif
        } // end of else if ( abandon )
      } // end of if( actAndReport || abandon || inform )
    } // end of if ( act || inform )
  }
  else
  {
    if (!fPushed)  { fNumberZeroSteps = 0; }
  }
  fLastMotherPhys = motherPhysical;

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

    if ( fValidExitNormal || fCalculatedExitNormal )
    {
      // Convention: fExitNormal is in the 'grand-mother' coordinate system
      fGrandMotherExitNormal = fExitNormal;
    }
    else
    {  
      // We must calculate the normal anyway (in order to have it if requested)
      //
      G4ThreeVector finalLocalPoint = fLastLocatedPointLocal
                                    + localDirection*Step;

      if (  fHistory.GetTopVolumeType() != kReplica )
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
          fChangedGrandMotherRefFrame = true;           
          fGrandMotherExitNormal = (*mRot).inverse() * exitNormalMotherFrame;
        }
        else
        {
          fGrandMotherExitNormal = exitNormalMotherFrame;
        }

        // Do not set fValidExitNormal -- this signifies
        // that the solid is convex!
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

    if ( fHistory.GetTopVolumeType() != kReplica )
      fCalculatedExitNormal = true;

    // Now transform it to the global reference frame !!
    //
    if( fValidExitNormal || fCalculatedExitNormal )
    {
      G4int depth = (G4int)fHistory.GetDepth();
      if( depth > 0 )
      {
        fExitNormalGlobalFrame = fHistory.GetTransform(depth-1)
                                .InverseTransformAxis( fGrandMotherExitNormal );
      }
      else
      {
        fExitNormalGlobalFrame = fGrandMotherExitNormal;
      }
    }
    else
    {
      fExitNormalGlobalFrame = G4ThreeVector( 0., 0., 0.);
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
    if( fVerbose > 5 )  { G4cout << G4endl; }
    if( Step == kInfinity )
    {
       G4cout << " Requested step= " << pCurrentProposedStepLength ;
       if( fVerbose > 5)  { G4cout << G4endl; }
    }
    G4cout << "  Safety = " << pNewSafety << G4endl;
  }
#endif

  fLastTriedStepComputation = true;

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
  fChangedGrandMotherRefFrame = false;
  fCalculatedExitNormal  = false;

  fExitNormal            = G4ThreeVector(0,0,0);
  fGrandMotherExitNormal = G4ThreeVector(0,0,0);
  fExitNormalGlobalFrame = G4ThreeVector(0,0,0);

  fPreviousSftOrigin     = G4ThreeVector(0,0,0);
  fPreviousSafety        = 0.0; 

  fNumberZeroSteps       = 0;

  fBlockedPhysicalVolume = nullptr;
  fBlockedReplicaNo      = -1;

  fLastLocatedPointLocal = G4ThreeVector( kInfinity, -kInfinity, 0.0 ); 
  fLocatedOutsideWorld   = false;

  fLastMotherPhys = nullptr;
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
  const G4int depth = (G4int)fHistory.GetDepth();
  for ( auto i = 1; i <= depth; ++i )
  {
    switch ( fHistory.GetVolumeType(i) )
    {
      case kNormal:
      case kExternal:
        break;
      case kReplica:
        freplicaNav.ComputeTransformation(fHistory.GetReplicaNo(i), fHistory.GetVolume(i));
        break;
      case kParameterised:
        G4VPhysicalVolume* current = fHistory.GetVolume(i);
        G4int replicaNo = fHistory.GetReplicaNo(i);
        G4VPVParameterisation* pParam = current->GetParameterisation();
        G4VSolid* pSolid = pParam->ComputeSolid(replicaNo, current);

        // Set up dimensions & transform in solid/physical volume
        //
        pSolid->ComputeDimensions(pParam, replicaNo, current);
        pParam->ComputeTransformation(replicaNo, current);

        G4TouchableHistory* pTouchable = nullptr;
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
        G4LogicalVolume* pLogical = current->GetLogicalVolume();
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
  G4VSolid* currentSolid = nullptr;
  G4LogicalVolume* candidateLogical;

  if ( fLastTriedStepComputation ) 
  {
    // use fLastLocatedPointLocal and next candidate volume
    //
    G4ThreeVector nextSolidExitNormal(0.,0.,0.);

    if( fEntering && (fBlockedPhysicalVolume!=0) ) 
    { 
      candidateLogical = fBlockedPhysicalVolume->GetLogicalVolume();
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
          G4ThreeVector daughterPointOwnLocal =
            MotherToDaughterTransform.TransformPoint( fLastStepEndPointLocal ); 

          // OK if it is a parameterised volume
          //
          EInside inSideIt; 
          G4bool onSurface;
          G4double safety = -1.0; 
          currentSolid = candidateLogical->GetSolid(); 
          inSideIt = currentSolid->Inside(daughterPointOwnLocal); 
          onSurface = (inSideIt == kSurface); 
          if( !onSurface ) 
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
            // First flip ( ExitNormal = -nextSolidExitNormal; )
            //  and then rotate the the normal to the frame of the mother (current volume)
            ExitNormal = MotherToDaughterTransform
                        .InverseTransformAxis( -nextSolidExitNormal );
            fCalculatedExitNormal = true;
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
      fCalculatedExitNormal = true;  // Should be true already
    }
    else  // i.e.  ( fBlockedPhysicalVolume == 0 )
    {
      *valid = false;
      G4Exception("G4Navigator::GetLocalExitNormal()",
                  "GeomNav0003", JustWarning, 
                  "Incorrect call to GetLocalSurfaceNormal." );
    }
  }
  else //  ( ! fLastTriedStepComputation ) i.e. last call was to Locate
  {
    if ( EnteredDaughterVolume() )
    {
      G4VSolid* daughterSolid = fHistory.GetTopVolume()->GetLogicalVolume()
                                                       ->GetSolid();
      ExitNormal = -(daughterSolid->SurfaceNormal(fLastLocatedPointLocal));
      if( std::fabs(ExitNormal.mag2()-1.0 ) > kToleranceNormalCheck )
      {
        G4ExceptionDescription desc;
        desc << " Parameters of solid: " << *daughterSolid
             << " Point for surface = " << fLastLocatedPointLocal << std::endl;
        G4Exception("G4Navigator::GetLocalExitNormal()",
                    "GeomNav0003", FatalException, desc,
                    "Surface Normal returned by Solid is not a Unit Vector." );
      }
      fCalculatedExitNormal = true;
      *valid = true;
    }
    else
    {
      if( fExitedMother )
      {
        ExitNormal = fGrandMotherExitNormal;
        *valid = true;
        fCalculatedExitNormal = true;
      }
      else  // We are not at a boundary. ExitNormal remains (0,0,0)
      { 
        *valid = false;
        fCalculatedExitNormal = false; 
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
G4Navigator::GetMotherToDaughterTransform( G4VPhysicalVolume* pEnteringPhysVol,
                                           G4int enteringReplicaNo,
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
    case kExternal:
      // Expect that nothing is needed to prepare the transformation.
      // It is stored already in the physical volume (placement)
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
                                 G4bool* pValid)
{
#ifdef G4DEBUG_NAVIGATION
  // Check Current point against expected 'local' value
  //
  if ( fLastTriedStepComputation ) 
  {
    G4ThreeVector ExpectedBoundaryPointLocal;

    const G4AffineTransform& GlobalToLocal = GetGlobalToLocalTransform(); 
    ExpectedBoundaryPointLocal =
      GlobalToLocal.TransformPoint( ExpectedBoundaryPointGlobal ); 

    // Add here:  Comparison against expected position,
    //            i.e. the endpoint of ComputeStep
  }
#endif
  
  return GetLocalExitNormal( pValid ); 
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
    if( std::fabs ( normMag2 - 1.0 ) < perThousand ) // was perMillion
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
       fVerbose = 4; 
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
       globalNormal = fHistory.GetTopTransform()
                     .InverseTransformAxis(localNormal);
    }
  }
  else
  {
    localNormal = GetLocalExitNormalAndCheck(IntersectPointGlobal,&validNormal);
    *pNormalCalculated = fCalculatedExitNormal;

#ifdef G4DEBUG_NAVIGATION
    usingStored = false;

    if( (!validNormal) && !fCalculatedExitNormal )
    {
      G4ExceptionDescription edN;
      edN << "  Calculated = " << fCalculatedExitNormal << G4endl;
      edN << "   Entering= "  << fEntering << G4endl;
      G4int oldVerbose = this->GetVerboseLevel();
      this->SetVerboseLevel(4);
      edN << "   State of Navigator: " << G4endl;
      edN << *this << G4endl;
      this->SetVerboseLevel( oldVerbose );
       
      G4Exception("G4Navigator::GetGlobalExitNormal()",
                  "GeomNav0003", JustWarning, edN,
                  "LocalExitNormalAndCheck() did not calculate Normal.");
     }
#endif
     
     G4double localMag2 = localNormal.mag2();
     if( validNormal && (std::fabs(localMag2-1.0)) > kToleranceNormalCheck )
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
     globalNormal = fHistory.GetTopTransform()
                   .InverseTransformAxis(localNormal);
  }

#ifdef G4DEBUG_NAVIGATION
  if( usingStored )
  {
    G4ThreeVector globalNormAgn; 

    localNormal = GetLocalExitNormalAndCheck(IntersectPointGlobal,&validNormal);
    
    globalNormAgn = fHistory.GetTopTransform()
                   .InverseTransformAxis(localNormal);
    
    // Check the value computed against fExitNormalGlobalFrame
    G4ThreeVector diffNorm = globalNormAgn - fExitNormalGlobalFrame;
    if( diffNorm.mag2() > kToleranceNormalCheck )
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

// ********************************************************************
// ComputeSafety
//
// It assumes that it will be 
//  i) called at the Point in the same volume as the EndPoint of the
//     ComputeStep.
// ii) after (or at the end of) ComputeStep OR after the relocation.
// ********************************************************************
//
G4double G4Navigator::ComputeSafety( const G4ThreeVector& pGlobalpoint,
                                     const G4double pMaxLength,
                                     const G4bool keepState)
{
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
  G4bool stayedOnEndpoint = distEndpointSq < sqr(kCarTolerance); 
  G4bool endpointOnSurface = fEnteredDaughter || fExitedMother;

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
      G4cout << "   ---- Exiting ComputeSafety  " << G4endl;
      PrintState();
      G4cout << "    Returned value of Safety is zero " << G4endl;
      G4cout.precision(oldcoutPrec);
    }
#endif
    return 0.0;
  }

  G4double newSafety = 0.0;

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
    G4VPhysicalVolume* motherPhysical = fHistory.GetTopVolume();
    G4LogicalVolume* motherLogical = motherPhysical->GetLogicalVolume();
    G4SmartVoxelHeader* pVoxelHeader = motherLogical->GetVoxelHeader();
    G4ThreeVector localPoint = ComputeLocalPoint(pGlobalpoint);

    if ( fHistory.GetTopVolumeType() != kReplica )
    {
      switch(CharacteriseDaughters(motherLogical))
      {
        case kNormal:
          if ( pVoxelHeader )
          {
            newSafety = fpVoxelSafety->ComputeSafety(localPoint,
                                             *motherPhysical, pMaxLength);
            // = VoxelNav().ComputeSafety(localPoint,fHistory,pMaxLength); // - Old method
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
        case kExternal:
          newSafety = fpExternalNav->ComputeSafety(localPoint, fHistory,
                                                   pMaxLength);
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
  G4long oldcoutPrec = G4cout.precision(4);
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
           << "  LastStepWasZero      = " <<  fLastStepWasZero    //  << G4endl
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
    if ( fBlockedPhysicalVolume == nullptr )
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

  G4ThreeVector OriginalGlobalpoint = fHistory.GetTopTransform().
                               InverseTransformPoint(fLastLocatedPointLocal); 

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
      G4long oldcoutPrec = G4cout.precision(8);
      G4long oldcerrPrec = G4cerr.precision(10);
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
// CheckOverlapsIterative
// ********************************************************************
//
G4bool G4Navigator::CheckOverlapsIterative(G4VPhysicalVolume* vol)
{
  // Check and report overlaps
  //
  G4bool foundOverlap = false;
  G4int  nPoints = 300000,  ntrials = 9, numOverlaps = 5;
  G4double  trialLength = 1.0 * CLHEP::centimeter;
  while ( ntrials-- > 0 && !foundOverlap )
  {
    if ( fVerbose > 1 )
    {
       G4cout << " ** Running overlap checks in volume "
              <<  vol->GetName()
              << " with length = " << trialLength << G4endl;
    }
    foundOverlap = vol->CheckOverlaps(nPoints, trialLength,
                                      fVerbose, numOverlaps);
    trialLength *= 0.1;
    if ( trialLength <= 1.0e-5 ) { numOverlaps= 1;}
  }
  return foundOverlap;
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

  G4long oldcoutPrec = os.precision(4);
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

#ifdef ALTERNATIVE_VOXEL_NAV
// ********************************************************************
// SetVoxelNavigation  -- alternative navigator for Voxel geom
// ********************************************************************
//
void G4Navigator::SetVoxelNavigation(G4VoxelNavigation* voxelNav)
{
  delete fpvoxelNav;
  fpvoxelNav = voxelNav;
}
#endif

// ********************************************************************
// InformLastStep:  Derived navigators can inform of its step
//                    - used to update fLastStepWasZero
// ********************************************************************
void  G4Navigator::InformLastStep(G4double lastStep, G4bool entersDaughtVol, G4bool exitsMotherVol )
{
  G4bool zeroStep = ( lastStep == 0.0 );
  fLocatedOnEdge   = fLastStepWasZero && zeroStep;  
  fLastStepWasZero = zeroStep;

  fExiting = exitsMotherVol;
  fEntering = entersDaughtVol;
}
