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
// $Id: G4PathFinder.cc 103219 2017-03-22 11:30:15Z gcosmo $
// GEANT4 tag $ Name:  $
// 
// class G4PathFinder Implementation
//
// Original author:  John Apostolakis,  April 2006
// 
// --------------------------------------------------------------------

#include <iomanip>

#include "G4PathFinder.hh"

#include "G4SystemOfUnits.hh"
#include "G4GeometryTolerance.hh"
#include "G4Navigator.hh"
#include "G4PropagatorInField.hh"
#include "G4TransportationManager.hh"
#include "G4MultiNavigator.hh"
#include "G4SafetyHelper.hh"

// Initialise the static instance of the singleton
//
G4ThreadLocal G4PathFinder* G4PathFinder::fpPathFinder=0;

// ----------------------------------------------------------------------------
// GetInstance()
//
// Retrieve the static instance of the singleton and create it if not existing
//
G4PathFinder* G4PathFinder::GetInstance()
{
   if( ! fpPathFinder )
   {
     fpPathFinder = new G4PathFinder; 
   }
   return fpPathFinder;
}

// ----------------------------------------------------------------------------
// GetInstanceIfExist()
//
// Retrieve the static instance pointer of the singleton
//
G4PathFinder* G4PathFinder::GetInstanceIfExist()
{
   return fpPathFinder;
}

// ----------------------------------------------------------------------------
// Constructor
//
G4PathFinder::G4PathFinder() 
  : fEndState( G4ThreeVector(), G4ThreeVector(), 0., 0., 0., 0., 0.),
    fFieldExertedForce(false),
    fRelocatedPoint(true),
    fLastStepNo(-1), fCurrentStepNo(-1),
    fVerboseLevel(0)
{
   fpMultiNavigator= new G4MultiNavigator(); 

   fpTransportManager= G4TransportationManager::GetTransportationManager();
   fpFieldPropagator = fpTransportManager->GetPropagatorInField();

   kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

   fNoActiveNavigators= 0; 
   G4ThreeVector  Big3Vector( kInfinity, kInfinity, kInfinity );
   fLastLocatedPosition= Big3Vector;
   fSafetyLocation= Big3Vector; 
   fPreSafetyLocation= Big3Vector;
   fPreStepLocation= Big3Vector;

   fPreSafetyMinValue=  -1.0; 
   fMinSafety_PreStepPt= -1.0; 
   fMinSafety_atSafLocation= -1.0; 
   fMinStep=   -1.0;
   fTrueMinStep= -1.0;
   fPreStepCenterRenewed= false;
   fNewTrack= false; 
   fNoGeometriesLimiting= 0; 

   for( G4int num=0; num< fMaxNav; ++num )
   {
      fpNavigator[num] =  0;   
      fLimitTruth[num] = false;
      fLimitedStep[num] = kUndefLimited;
      fCurrentStepSize[num] = -1.0; 
      fLocatedVolume[num] = 0;
      fPreSafetyValues[num]= -1.0; 
      fCurrentPreStepSafety[num] = -1.0;
      fNewSafetyComputed[num]= -1.0; 
   }
}

// ----------------------------------------------------------------------------
// Destructor
//
G4PathFinder::~G4PathFinder() 
{
   delete fpMultiNavigator;
   fpPathFinder = 0;
}

// ----------------------------------------------------------------------------
//
void
G4PathFinder::EnableParallelNavigation(G4bool enableChoice)
{
   G4Navigator *navigatorForPropagation=0, *massNavigator=0;

   massNavigator= fpTransportManager->GetNavigatorForTracking(); 
   if( enableChoice )
   {
      navigatorForPropagation= fpMultiNavigator;

      // Enable SafetyHelper to use PF
      //
      fpTransportManager->GetSafetyHelper()->EnableParallelNavigation(true);
   }
   else
   {
      navigatorForPropagation= massNavigator;
       
      // Disable SafetyHelper to use PF
      //
      fpTransportManager->GetSafetyHelper()->EnableParallelNavigation(false);
   }
   fpFieldPropagator->SetNavigatorForPropagating(navigatorForPropagation);
}

// ----------------------------------------------------------------------------
//
G4double 
G4PathFinder::ComputeStep( const G4FieldTrack &InitialFieldTrack, 
                                 G4double     proposedStepLength,
                                 G4int        navigatorNo, 
                                 G4int        stepNo,       // find next step 
                                 G4double     &pNewSafety,  // for this geom 
                                 ELimited     &limitedStep, 
                                 G4FieldTrack &EndState,
                                 G4VPhysicalVolume* currentVolume)
{
  G4double  possibleStep= -1.0; 

#ifdef G4DEBUG_PATHFINDER
  if( fVerboseLevel > 2 )
  { 
    G4cout << " -------------------------" <<  G4endl;
    G4cout << " G4PathFinder::ComputeStep - entered " << G4endl;
    G4cout << "   - stepNo = "  << std::setw(4) << stepNo  << " "
           << " navigatorId = " << std::setw(2) << navigatorNo  << " "
           << " proposed step len = " << proposedStepLength << " " << G4endl;
    G4cout << " PF::ComputeStep requested step " 
           << " from " << InitialFieldTrack.GetPosition()
           << " dir  " << InitialFieldTrack.GetMomentumDirection() << G4endl;
  }
#endif
#ifdef G4VERBOSE
  if( navigatorNo >= fNoActiveNavigators )
  {
    std::ostringstream message;
    message << "Bad Navigator ID !" << G4endl
            << "        Requested Navigator ID = " << navigatorNo << G4endl
            << "        Number of active navigators = " << fNoActiveNavigators;
    G4Exception("G4PathFinder::ComputeStep()", "GeomNav0002",
                FatalException, message); 
  }
#endif

  if( fNewTrack || (stepNo != fLastStepNo)  )
  {
    // This is a new track or a new step, so we must make the step
    // ( else we can simply retrieve its results for this Navigator Id )    

    G4FieldTrack currentState= InitialFieldTrack;

    fCurrentStepNo = stepNo; 

    // Check whether a process shifted the position 
    // since the last step -- by physics processes
    //
    G4ThreeVector newPosition = InitialFieldTrack.GetPosition();   
    G4ThreeVector moveVector= newPosition - fLastLocatedPosition; 
    G4double moveLenSq= moveVector.mag2(); 
    if( moveLenSq > kCarTolerance * kCarTolerance )
    { 
       G4ThreeVector newDirection = InitialFieldTrack.GetMomentumDirection();   
#ifdef G4DEBUG_PATHFINDER
       if( fVerboseLevel > 2 )
       { 
          G4double moveLen= std::sqrt( moveLenSq ); 
          G4cout << " G4PathFinder::ComputeStep : Point moved since last step " 
                 << " -- at step # = " << stepNo << G4endl
                 << " by " << moveLen  << " to " << newPosition << G4endl;      
       } 
#endif
       MovePoint();  // Unintentional changed -- ????

       // Relocate to cope with this move -- else could abort !?
       //
       Locate( newPosition, newDirection ); 
    }

    // Check whether the particle have an (EM) field force exerting upon it
    //
    G4double particleCharge=  currentState.GetCharge(); 

    G4FieldManager* fieldMgr=0;
    G4bool          fieldExertsForce = false ;
    if( (particleCharge != 0.0) )
    {
        fieldMgr= fpFieldPropagator->FindAndSetFieldManager( currentVolume );

        // Protect for case where field manager has no field (= field is zero)
        //
        fieldExertsForce = (fieldMgr != 0) 
                        && (fieldMgr->GetDetectorField() != 0);
    }
    fFieldExertedForce = fieldExertsForce;  // Store for use in later calls
                                            // referring to this 'step'.

    fNoGeometriesLimiting= -1;  // At start of track, no process limited step
    if( fieldExertsForce )
    {
       DoNextCurvedStep( currentState, proposedStepLength, currentVolume ); 
       //--------------
    }else{
       DoNextLinearStep( currentState, proposedStepLength ); 
       //--------------
    }
    fLastStepNo= stepNo; 

#ifdef  G4DEBUG_PATHFINDER
    if ( (fNoGeometriesLimiting < 0)
      || (fNoGeometriesLimiting > fNoActiveNavigators) )
    {
      std::ostringstream message;
      message << "Number of geometries limiting the step not set." << G4endl
              << "        Number of geometries limiting step = "
              << fNoGeometriesLimiting;
      G4Exception("G4PathFinder::ComputeStep()", 
                  "GeomNav0002", FatalException, message); 
    }
#endif
  }
#ifdef G4DEBUG_PATHFINDER      
  else
  {
     const G4double checkTolerance = 1.0e-9; 
     if( proposedStepLength < fTrueMinStep * ( 1.0 + checkTolerance) )  // For 2nd+ geometry 
     { 
       std::ostringstream message;
       message.precision( 12 ); 
       message << "Problem in step size request." << G4endl
               << "        Being requested to make a step which is shorter"
               << " than the minimum Step " << G4endl
               << "        already computed for any Navigator/geometry during"
               << " this tracking-step: " << G4endl
               << "        This could happen due to an error in process ordering."
               << G4endl
               << "        Check that all physics processes are registered"
               << "        before all processes with a navigator/geometry."
               << G4endl
               << "        If using pre-packaged physics list and/or"
               << "        functionality, please report this error."
               << G4endl << G4endl
               << "        Additional information for problem: "  << G4endl
               << "        Steps request/proposed = " << proposedStepLength
               << G4endl
               << "        MinimumStep (true) = " << fTrueMinStep
               << G4endl
               << "        MinimumStep (navraw)  = " << fMinStep
               << G4endl
               << "        Navigator raw return value" << G4endl
               << "        Requested step now = " << proposedStepLength
               << G4endl
               << "        Difference min-req (absolute) = "
               << fTrueMinStep-proposedStepLength << G4endl
               << "        Relative (to max of two) = " 
               << (fTrueMinStep-proposedStepLength)
                  / std::max(proposedStepLength, fTrueMinStep) << G4endl
               << "     -- Step info> stepNo= " << stepNo
               << " last= " << fLastStepNo 
               << " newTr= " << fNewTrack << G4endl;
        G4Exception("G4PathFinder::ComputeStep()", 
                    "GeomNav0003", FatalException, message);
     }
     else
     { 
        // This is neither a new track nor a new step -- just another 
        // client accessing information for the current track, step 
        // We will simply retrieve the results of the synchronous
        // stepping for this Navigator Id below.
        //
        if( fVerboseLevel > 1 )
        { 
           G4cout << " G4P::CS -> Not calling DoNextLinearStep: " 
                  << " stepNo= " << stepNo << " last= " << fLastStepNo 
                  << " new= " << fNewTrack << " Step already done" << G4endl; 
        }
     } 
  }
#endif

  fNewTrack= false; 

  // Prepare the information to return

  pNewSafety  = fCurrentPreStepSafety[ navigatorNo ]; 
  limitedStep = fLimitedStep[ navigatorNo ];
  fRelocatedPoint= false;

  possibleStep= std::min(proposedStepLength, fCurrentStepSize[ navigatorNo ]);
  EndState = fEndState;  //  now corrected for smaller step, if needed

#ifdef G4DEBUG_PATHFINDER
  if( fVerboseLevel > 0 )
  { 
    G4cout << " G4PathFinder::ComputeStep returns "
           << fCurrentStepSize[ navigatorNo ]
           << " for Navigator " << navigatorNo 
           << " Limited step = " << limitedStep 
           << " Safety(mm) = " << pNewSafety / mm 
           << G4endl; 
  }
#endif

  return possibleStep;
}

// ----------------------------------------------------------------------

void
G4PathFinder::PrepareNewTrack( const G4ThreeVector& position, 
                               const G4ThreeVector& direction,
                               G4VPhysicalVolume*  massStartVol)
{
  // Key purposes:
  //   - Check and cache set of active navigators
  //   - Reset state for new track

  G4int num=0; 

  EnableParallelNavigation(true); 
    // Switch PropagatorInField to use MultiNavigator

  fpTransportManager->GetSafetyHelper()->InitialiseHelper(); 
    // Reinitialise state of safety helper -- avoid problems with overlaps

  fNewTrack= true; 
  this->MovePoint();   // Signal further that the last status is wiped

  fpFieldPropagator->PrepareNewTrack();   // Inform field propagator of new track
  
  // Message the G4NavigatorPanel / Dispatcher to find active navigators
  //
  std::vector<G4Navigator*>::iterator pNavigatorIter; 

  fNoActiveNavigators=  fpTransportManager-> GetNoActiveNavigators();
  if( fNoActiveNavigators > fMaxNav )
  {
    std::ostringstream message;
    message << "Too many active Navigators / worlds." << G4endl
            << "        Transportation Manager has "
            << fNoActiveNavigators << " active navigators." << G4endl
            << "        This is more than the number allowed = "
            << fMaxNav << " !";
    G4Exception("G4PathFinder::PrepareNewTrack()", "GeomNav0002",  
                FatalException, message); 
  }

  fpMultiNavigator->PrepareNavigators(); 
  //------------------------------------

  pNavigatorIter= fpTransportManager->GetActiveNavigatorsIterator();
  for( num=0; num< fNoActiveNavigators; ++pNavigatorIter,++num )
  {
     // Keep information in C-array ... for creating touchables - at least

     fpNavigator[num] =  *pNavigatorIter;   
     fLimitTruth[num] = false;
     fLimitedStep[num] = kDoNot;
     fCurrentStepSize[num] = 0.0; 
     fLocatedVolume[num] = 0; 
  }
  fNoGeometriesLimiting= 0;  // At start of track, no process limited step

  // In case of one geometry, the tracking will have done the locating!!

  if( fNoActiveNavigators > 1 )
  {
     Locate( position, direction, false );   
  }
  else
  {
     // Update state -- depending on the tracking's call to Mass Navigator

     fLastLocatedPosition= position; 
     fLocatedVolume[0]= massStartVol; // This information must be given
                                      // by transportation
     fLimitedStep[0]   = kDoNot; 
     fCurrentStepSize[0] = 0.0;
  }

  // Reset Safety Information -- as in case of overlaps this can cause
  // inconsistencies ...
  //
  fMinSafety_PreStepPt= fPreSafetyMinValue= fMinSafety_atSafLocation= 0.0; 
 
  for( num=0; num< fNoActiveNavigators; ++num )
  {
     fPreSafetyValues[num]= 0.0; 
     fNewSafetyComputed[num]= 0.0; 
     fCurrentPreStepSafety[num] = 0.0;
  }

  // The first location for each Navigator must be non-relative
  // or else call ResetStackAndState() for each Navigator

  fRelocatedPoint= false; 
}


void G4PathFinder::EndTrack()
     // Signal end of tracking of current track.  
     //   Reset TransportationManager to use 'ordinary' Navigator
     //   Reset internal state, if needed
{
  EnableParallelNavigation(false);  // Else it will be continue to be used
}

void G4PathFinder::ReportMove( const G4ThreeVector& OldVector, 
                               const G4ThreeVector& NewVector, 
                               const G4String& Quantity ) const
{
    G4ThreeVector moveVec = ( NewVector - OldVector );

    G4int prc= G4cerr.precision(12); 
    std::ostringstream message;
    message << "Endpoint moved between value returned by ComputeStep()"
            << " and call to Locate(). " << G4endl
            << "          Change of " << Quantity << " is "
            << moveVec.mag() / mm << " mm long" << G4endl
            << "          and its vector is "
            << (1.0/mm) * moveVec << " mm " << G4endl
            << "          Endpoint of ComputeStep() was " << OldVector << G4endl
            << "          and current position to locate is " << NewVector;
    G4Exception("G4PathFinder::ReportMove()", "GeomNav1002",  
                JustWarning, message); 
    G4cerr.precision(prc); 
}

void
G4PathFinder::Locate( const   G4ThreeVector& position, 
                      const   G4ThreeVector& direction,
                      G4bool  relative)
{
  // Locate the point in each geometry

  std::vector<G4Navigator*>::iterator pNavIter=
     fpTransportManager->GetActiveNavigatorsIterator(); 

  G4ThreeVector lastEndPosition= fEndState.GetPosition(); 
  G4ThreeVector moveVec = (position - lastEndPosition );
  G4double      moveLenSq= moveVec.mag2();
  if( (!fNewTrack) && (!fRelocatedPoint)
   && ( moveLenSq> 10*kCarTolerance*kCarTolerance ) )
  {
     ReportMove( position, lastEndPosition, "Position" ); 
  }
  fLastLocatedPosition= position; 

#ifdef G4DEBUG_PATHFINDER
  if( fVerboseLevel > 2 )
  {
    G4cout << G4endl; 
    G4cout << " G4PathFinder::Locate : entered " << G4endl;
    G4cout << " --------------------   -------" <<  G4endl;
    G4cout << "   Locating at position " << position
           << "  with direction " << direction 
           << "  relative= " << relative << G4endl;
    if ( (fVerboseLevel > 1) || ( moveLenSq > 0.0) )
    { 
       G4cout << "  lastEndPosition = " << lastEndPosition
              << "  moveVec = " << moveVec
              << "  newTr = " << fNewTrack 
              << "  relocated = " << fRelocatedPoint << G4endl;
    }

    G4cout << " Located at " << position ; 
    if( fNoActiveNavigators > 1 )  { G4cout << G4endl; }
  }
#endif

  for ( G4int num=0; num< fNoActiveNavigators ; ++pNavIter,++num )
  {
     //  ... who limited the step ....

     if( fLimitTruth[num] ) { (*pNavIter)->SetGeometricallyLimitedStep(); }

     G4VPhysicalVolume *pLocated= 
     (*pNavIter)->LocateGlobalPointAndSetup( position, &direction,
                                             relative,  
                                             false);   
     // Set the state related to the location
     //
     fLocatedVolume[num] = pLocated; 

     // Clear state related to the step
     //
     fLimitedStep[num]   = kDoNot; 
     fCurrentStepSize[num] = 0.0;      
    
#ifdef G4DEBUG_PATHFINDER
     if( fVerboseLevel > 2 )
     {
       G4cout << " - In world " << num << " geomLimStep= " << fLimitTruth[num]
              << "  gives volume= " << pLocated ; 
       if( pLocated )
       { 
         G4cout << "  name = '" << pLocated->GetName() << "'"; 
         G4cout << " - CopyNo= " << pLocated->GetCopyNo(); 
       } 
       G4cout  << G4endl; 
     }
#endif
  }

  fRelocatedPoint= false;
}

void G4PathFinder::ReLocate( const G4ThreeVector& position )
{
  // Locate the point in each geometry

  std::vector<G4Navigator*>::iterator pNavIter=
    fpTransportManager->GetActiveNavigatorsIterator(); 

#ifdef G4DEBUG_PATHFINDER

  // Check that this relocation does not violate safety
  //   - at endpoint (computed from start point) AND
  //   - at last safety location  (likely just called)

  G4ThreeVector lastEndPosition= fEndState.GetPosition();

  // Calculate end-point safety ...
  //
  G4double      DistanceStartEnd= (lastEndPosition - fPreStepLocation).mag();
  G4double      endPointSafety_raw = fMinSafety_PreStepPt - DistanceStartEnd; 
  G4double      endPointSafety_Est1 = std::max( 0.0, endPointSafety_raw ); 

  // ... and check move from endpoint against this endpoint safety
  //
  G4ThreeVector moveVecEndPos  = position - lastEndPosition;
  G4double      moveLenEndPosSq = moveVecEndPos.mag2(); 

  // Check that move from endpoint of last step is within safety
  // -- or check against last location or relocation ?? 
  //
  G4ThreeVector moveVecSafety=  position - fSafetyLocation; 
  G4double      moveLenSafSq=   moveVecSafety.mag2();

  G4double distCheckEnd_sq= ( moveLenEndPosSq - endPointSafety_Est1 
                                               *endPointSafety_Est1 ); 
  G4double distCheckSaf_sq=   ( moveLenSafSq -  fMinSafety_atSafLocation
                                               *fMinSafety_atSafLocation ); 

  G4bool longMoveEnd = distCheckEnd_sq > 0.0; 
  G4bool longMoveSaf = distCheckSaf_sq > 0.0; 

  G4double revisedSafety= 0.0;

  if( (!fNewTrack) && ( longMoveEnd && longMoveSaf ) )
  {  
     // Recompute ComputeSafety for end position
     //
     revisedSafety= ComputeSafety(lastEndPosition); 

     const G4double kRadTolerance =
       G4GeometryTolerance::GetInstance()->GetRadialTolerance();
     const G4double cErrorTolerance=1e-12;   
       // Maximum relative error from roundoff of arithmetic 

     G4double  distCheckRevisedEnd= moveLenEndPosSq-revisedSafety*revisedSafety;

     G4bool  longMoveRevisedEnd=  ( distCheckRevisedEnd > 0. ) ; 

     G4double  moveMinusSafety= 0.0; 
     G4double  moveLenEndPosition= std::sqrt( moveLenEndPosSq );
     moveMinusSafety = moveLenEndPosition - revisedSafety; 

     if ( longMoveRevisedEnd && (moveMinusSafety > 0.0 )
       && ( revisedSafety > 0.0 ) )
     {
        // Take into account possibility of roundoff error causing
        // this apparent move further than safety

        if( fVerboseLevel > 0 )
        {
           G4cout << " G4PF:Relocate> Ratio to revised safety is " 
                  << std::fabs(moveMinusSafety)/revisedSafety << G4endl;
        }

        G4double  absMoveMinusSafety= std::fabs(moveMinusSafety);
        G4bool smallRatio= absMoveMinusSafety < kRadTolerance * revisedSafety ; 
        G4double maxCoordPos = std::max( 
                                      std::max( std::fabs(position.x()), 
                                                std::fabs(position.y())), 
                                      std::fabs(position.z()) );
        G4bool smallValue= absMoveMinusSafety < cErrorTolerance * maxCoordPos;
        if( ! (smallRatio || smallValue) )
        {
           G4cout << " G4PF:Relocate> Ratio to revised safety is " 
                  << std::fabs(moveMinusSafety)/revisedSafety << G4endl;
           G4cout << " Difference of move and safety is not very small."
                  << G4endl;
        }
        else
        {
          moveMinusSafety = 0.0; 
          longMoveRevisedEnd = false;   // Numerical issue -- not too long!

          G4cout << " Difference of move & safety is very small in magnitude, "
                 << absMoveMinusSafety << G4endl;
          if( smallRatio )
          {
            G4cout << " ratio to safety " << revisedSafety 
                   << " is " <<  absMoveMinusSafety / revisedSafety
                   << "smaller than " << kRadTolerance << " of safety ";
          }
          else
          {
            G4cout << " as fraction " << absMoveMinusSafety / maxCoordPos 
                   << " of position vector max-coord " << maxCoordPos
                   << " smaller than " << cErrorTolerance ;
          }
          G4cout << " -- reset moveMinusSafety to "
                 << moveMinusSafety << G4endl;
        }
     }

     if ( longMoveEnd && longMoveSaf
       && longMoveRevisedEnd && (moveMinusSafety>0.0) )
     { 
        G4int oldPrec= G4cout.precision(9); 
        std::ostringstream message;
        message << "ReLocation is further than end-safety value." << G4endl
                << " Moved from last endpoint by " << moveLenEndPosition 
                << " compared to end safety (from preStep point) = " 
                << endPointSafety_Est1 << G4endl
                << "  --> last PreSafety Location was " << fPreSafetyLocation
                << G4endl
                << "       safety value =  " << fPreSafetyMinValue << G4endl
                << "  --> last PreStep Location was " << fPreStepLocation
                << G4endl
                << "       safety value =  " << fMinSafety_PreStepPt << G4endl
                << "  --> last EndStep Location was " << lastEndPosition
                << G4endl
                << "       safety value =  " << endPointSafety_Est1 
                << " raw-value = " << endPointSafety_raw << G4endl
                << "  --> Calling again at this endpoint, we get "
                <<  revisedSafety << " as safety value."  << G4endl
                << "  --> last position for safety " << fSafetyLocation
                << G4endl
                << "       its safety value =  " << fMinSafety_atSafLocation
                << G4endl
                << "       move from safety location = "
                << std::sqrt(moveLenSafSq) << G4endl
                << "         again= " << moveVecSafety.mag() << G4endl
                << "       safety - Move-from-end= " 
                << revisedSafety - moveLenEndPosition
                << " (negative is Bad.)" << G4endl
                << " Debug:  distCheckRevisedEnd = "
                << distCheckRevisedEnd;
        ReportMove( lastEndPosition, position, "Position" ); 
        G4Exception("G4PathFinder::ReLocate", "GeomNav0003", 
                    FatalException, message); 
        G4cout.precision(oldPrec); 
    }
  }

  if( fVerboseLevel > 2 )
  {
    G4cout << G4endl; 
    G4cout << " G4PathFinder::ReLocate : entered " << G4endl;
    G4cout << " ----------------------   -------" <<  G4endl;
    G4cout << "  *Re*Locating at position " << position  << G4endl; 
      // << "  with direction " << direction 
      // << "  relative= " << relative << G4endl;
    if ( (fVerboseLevel > -1) || ( moveLenEndPosSq > 0.0) )
    {
       G4cout << "  lastEndPosition = " << lastEndPosition
              << "  moveVec from step-end = " << moveVecEndPos
              << "  is new Track = " << fNewTrack 
              << "  relocated = " << fRelocatedPoint << G4endl;
    }
  }
#endif // G4DEBUG_PATHFINDER

  for ( G4int num=0; num< fNoActiveNavigators ; ++pNavIter,++num )
  {
     //  ... none limited the step

     (*pNavIter)->LocateGlobalPointWithinVolume( position ); 

     // Clear state related to the step
     //
     fLimitedStep[num]   = kDoNot; 
     fCurrentStepSize[num] = 0.0;      
     fLimitTruth[num] = false;   
  }

  fLastLocatedPosition= position; 
  fRelocatedPoint= false;

#ifdef G4DEBUG_PATHFINDER
  if( fVerboseLevel > 2 )
  {
    G4cout << " G4PathFinder::ReLocate : exiting " 
           << "  at position " << fLastLocatedPosition << G4endl << G4endl;
  }
#endif
}

// -----------------------------------------------------------------------------

G4double  G4PathFinder::ComputeSafety( const G4ThreeVector& position )
{
    // Recompute safety for the relevant point

   G4double minSafety= kInfinity; 
  
   std::vector<G4Navigator*>::iterator pNavigatorIter;
   pNavigatorIter= fpTransportManager->GetActiveNavigatorsIterator();

   for( G4int num=0; num<fNoActiveNavigators; ++pNavigatorIter,++num )
   {
      G4double safety = (*pNavigatorIter)->ComputeSafety( position, DBL_MAX, true );
      if( safety < minSafety ) { minSafety = safety; } 
      fNewSafetyComputed[num]= safety;
   } 

   fSafetyLocation= position;
   fMinSafety_atSafLocation = minSafety;

#ifdef G4DEBUG_PATHFINDER
   if( fVerboseLevel > 1 )
   { 
     G4cout << " G4PathFinder::ComputeSafety - returns " 
            << minSafety << " at location " << position << G4endl;
   }
#endif
   return minSafety; 
}


// -----------------------------------------------------------------------------

G4TouchableHandle 
G4PathFinder::CreateTouchableHandle( G4int navId ) const
{
#ifdef G4DEBUG_PATHFINDER
  if( fVerboseLevel > 2 )
  {
    G4cout << "G4PathFinder::CreateTouchableHandle : navId = "
           << navId << " -- " << GetNavigator(navId) << G4endl;
  }
#endif

  G4TouchableHistory* touchHist;
  touchHist= GetNavigator(navId) -> CreateTouchableHistory(); 

  G4VPhysicalVolume* locatedVolume= fLocatedVolume[navId]; 
  if( locatedVolume == 0 )
  {
     // Workaround to ensure that the touchable is fixed !! // TODO: fix

     touchHist->UpdateYourself( locatedVolume, touchHist->GetHistory() );
  }
 
#ifdef G4DEBUG_PATHFINDER
  if( fVerboseLevel > 2 )
  {   
    G4String VolumeName("None"); 
    if( locatedVolume ) { VolumeName= locatedVolume->GetName(); }
    G4cout << " Touchable History created at address " << touchHist
           << "  volume = " << locatedVolume << " name= " << VolumeName
           << G4endl;
  }
#endif

  return G4TouchableHandle(touchHist); 
}

G4double
G4PathFinder::DoNextLinearStep( const G4FieldTrack &initialState,
                                      G4double      proposedStepLength )
{
  std::vector<G4Navigator*>::iterator pNavigatorIter;
  G4double safety= 0.0, step=0.0;
  G4double minSafety= kInfinity, minStep;

  const G4int IdTransport= 0;  // Id of Mass Navigator !!
  G4int num=0; 

#ifdef G4DEBUG_PATHFINDER
  if( fVerboseLevel > 2 )
  {
    G4cout << " G4PathFinder::DoNextLinearStep : entered " << G4endl;
    G4cout << "   Input field track= " << initialState << G4endl;
    G4cout << "   Requested step= " << proposedStepLength << G4endl;
  }
#endif

  G4ThreeVector initialPosition= initialState.GetPosition(); 
  G4ThreeVector initialDirection= initialState.GetMomentumDirection();
  
  G4ThreeVector OriginShift = initialPosition - fPreSafetyLocation;
  G4double      MagSqShift  = OriginShift.mag2() ;
  G4double      MagShift;  // Only given value if it larger than minimum safety

  // Potential optimisation using Maximum Value of safety!
  // if( MagSqShift >= sqr(fPreSafetyMaxValue ) ){ 
  //   MagShift= kInfinity;   // Not a useful value -- all will not use/ignore
  // else
  //  MagShift= std::sqrt(MagSqShift) ;

  MagShift= std::sqrt(MagSqShift) ;

#ifdef G4PATHFINDER_OPTIMISATION

  G4double fullSafety;  // For all geometries, for prestep point

  if( MagSqShift >= sqr(fPreSafetyMinValue ) )
  {
     fullSafety = 0.0 ;     
  }
  else
  {
     fullSafety = fPreSafetyMinValue - MagShift;
  }
  if( proposedStepLength < fullSafety ) 
  {
     // Move is smaller than all safeties
     //  -> so we do not have to move the safety center

     fPreStepCenterRenewed= false;

     for( num=0; num< fNoActiveNavigators; ++num )
     {
        fCurrentStepSize[num]= kInfinity; 
        safety = std::max( 0.0,  fPreSafetyValues[num] - MagShift); 
        minSafety= std::min( safety, minSafety ); 
        fCurrentPreStepSafety[num]= safety; 
     }
     minStep= kInfinity;

#ifdef G4DEBUG_PATHFINDER
     if( fVerboseLevel > 2 )
     {
       G4cout << "G4PathFinder::DoNextLinearStep : Quick Stepping. " << G4endl
               << " proposedStepLength " <<  proposedStepLength
               << " < (full) safety = " << fullSafety 
               << " at " << initialPosition 
               << G4endl;
     }
#endif
  }
  else
#endif   // End of G4PATHFINDER_OPTIMISATION 1
  {
     // Move is larger than at least one of the safeties
     //  -> so we must move the safety center!

     fPreStepCenterRenewed= true;
     pNavigatorIter= fpTransportManager-> GetActiveNavigatorsIterator();

     minStep= kInfinity;  // Not proposedStepLength; 

     for( num=0; num< fNoActiveNavigators; ++pNavigatorIter,++num ) 
     {
        safety = std::max( 0.0,  fPreSafetyValues[num] - MagShift); 

#ifdef G4PATHFINDER_OPTIMISATION
        if( proposedStepLength <= safety )  // Should be just < safety ?
        {
           // The Step is guaranteed to be taken

           step= kInfinity;    //  ComputeStep Would return this

#ifdef G4DEBUG_PATHFINDER
           G4cout.precision(8); 
           G4cout << "PathFinder::ComputeStep> small proposed step = "
                  << proposedStepLength
                  << " <=  safety = " << safety << " for nav " << num 
                  << " Step fully taken. " << G4endl;
#endif
        }
        else
#endif   // End of G4PATHFINDER_OPTIMISATION 2
        {
#ifdef G4DEBUG_PATHFINDER
           G4double previousSafety= safety; 
#endif
           step= (*pNavigatorIter)->ComputeStep( initialPosition, 
                                                 initialDirection,
                                                 proposedStepLength,
                                                 safety ); 
           minStep  = std::min( step,  minStep);

           //  TODO: consider whether/how to reduce the proposed step 
           //        to the latest minStep value - to reduce calculations

#ifdef G4DEBUG_PATHFINDER
           if( fVerboseLevel > 0)
           {
             G4cout.precision(8); 
             G4cout << "PathFinder::ComputeStep> long  proposed step = "
                    << proposedStepLength
                    << "  >  safety = " << previousSafety
                    << " for nav " << num 
                    << " .  New safety = " << safety << " step= " << step
                    << G4endl;      
           } 
#endif
        }
        fCurrentStepSize[num] = step; 

        // Save safety value, must be done for all geometries "together"
        // (even if not recomputed using call to ComputeStep)
        // since they share the fPreSafetyLocation

        fPreSafetyValues[num]= safety; 
        fCurrentPreStepSafety[num]= safety; 

        minSafety= std::min( safety, minSafety ); 
           
#ifdef G4DEBUG_PATHFINDER
        if( fVerboseLevel > 2 )
        {
          G4cout << "G4PathFinder::DoNextLinearStep : Navigator ["
                 << num << "] -- step size " << step << G4endl;
        }
#endif
     }

     // Only change these when safety is recalculated
     // it is good/relevant only for safety calculations

     fPreSafetyLocation=  initialPosition; 
     fPreSafetyMinValue=  minSafety;
  } // end of else for  if( proposedStepLength <= fullSafety)

  // For use in Relocation, need PreStep point location, min-safety
  //
  fPreStepLocation= initialPosition; 
  fMinSafety_PreStepPt= minSafety; 

  fMinStep=   minStep; 

  if( fMinStep == kInfinity )
  {
     minStep = proposedStepLength;   //  Use this below for endpoint !!
  }
  fTrueMinStep = minStep;

  // Set the EndState

  G4ThreeVector endPosition;

  fEndState= initialState; 
  endPosition= initialPosition + minStep * initialDirection ; 

#ifdef G4DEBUG_PATHFINDER
  if( fVerboseLevel > 1 )
  {
    G4cout << "G4PathFinder::DoNextLinearStep : "
           << " initialPosition = " << initialPosition 
           << " and endPosition = " << endPosition<< G4endl;
  }
#endif

  fEndState.SetPosition( endPosition ); 
  fEndState.SetProperTimeOfFlight( -1.000 );   // Not defined YET

  if( fNoActiveNavigators == 1 )
  { 
     G4bool transportLimited = (fMinStep!= kInfinity); 
     fLimitTruth[IdTransport] = transportLimited; 
     fLimitedStep[IdTransport] = transportLimited ? kUnique : kDoNot;

     // Set fNoGeometriesLimiting - as WhichLimited does
     fNoGeometriesLimiting = transportLimited ? 1 : 0;  
  }
  else
  {
     WhichLimited(); 
  }

#ifdef G4DEBUG_PATHFINDER
  if( fVerboseLevel > 2 )
  {
    G4cout << " G4PathFinder::DoNextLinearStep : exits returning "
           << minStep << G4endl;
    G4cout << "   Endpoint values = " << fEndState << G4endl;
    G4cout << G4endl;
  }
#endif

  return minStep;
}

void G4PathFinder::WhichLimited()
{
  // Flag which processes limited the step

  G4int num=-1, last=-1; 
  G4int noLimited=0; 
  ELimited shared= kSharedOther; 

  const G4int IdTransport= 0;  // Id of Mass Navigator !!

  // Assume that [IdTransport] is Mass / Transport
  //
  G4bool transportLimited = (fCurrentStepSize[IdTransport] == fMinStep)
                           && ( fMinStep!= kInfinity) ; 

  if( transportLimited )  { 
     shared= kSharedTransport;
  }

  for ( num= 0; num < fNoActiveNavigators; num++ )
  { 
    G4bool limitedStep;

    G4double step= fCurrentStepSize[num]; 

    limitedStep = ( std::fabs(step - fMinStep) < kCarTolerance ) 
                 && ( step != kInfinity); 

    fLimitTruth[ num ] = limitedStep; 
    if( limitedStep )
    {
      noLimited++;  
      fLimitedStep[num] = shared;
      last= num; 
    }
    else
    {
      fLimitedStep[num] = kDoNot;
    }
  }
  fNoGeometriesLimiting= noLimited;  // Save # processes limiting step

  if( (last > -1) && (noLimited == 1 ) )
  {
    fLimitedStep[ last ] = kUnique; 
  }

#ifdef G4DEBUG_PATHFINDER
  if( fVerboseLevel > 1 )
  {
    PrintLimited();   // --> for tracing 
    if( fVerboseLevel > 4 ) {
      G4cout << " G4PathFinder::WhichLimited - exiting. " << G4endl;
    }
  }
#endif
}

void G4PathFinder::PrintLimited()
{
  // Report results -- for checking   

  G4cout << "G4PathFinder::PrintLimited reports: " ; 
  G4cout << "  Minimum step (true)= " << fTrueMinStep 
         << "  reported min = " << fMinStep 
         << G4endl; 
  if(  (fCurrentStepNo <= 2) || (fVerboseLevel>=2) )
  {
    G4cout << std::setw(5) << " Step#"  << " "
           << std::setw(5) << " NavId"  << " "
           << std::setw(12) << " step-size " << " "
           << std::setw(12) << " raw-size "  << " "
           << std::setw(12) << " pre-safety " << " " 
           << std::setw(15) << " Limited / flag"  << " "
           << std::setw(15) << "  World "  << " "
           << G4endl;  
  }
  G4int num;
  for ( num= 0; num < fNoActiveNavigators; num++ )
  { 
    G4double rawStep = fCurrentStepSize[num]; 
    G4double stepLen = fCurrentStepSize[num]; 
    if( stepLen > fTrueMinStep )
    { 
      stepLen = fTrueMinStep;     // did not limit (went as far as asked)
    }
    G4int oldPrec= G4cout.precision(9); 

    G4cout << std::setw(5) << fCurrentStepNo  << " " 
           << std::setw(5) << num  << " "
           << std::setw(12) << stepLen << " "
           << std::setw(12) << rawStep << " "
           << std::setw(12) << fCurrentPreStepSafety[num] << " "
           << std::setw(5) << (fLimitTruth[num] ? "YES" : " NO") << " ";
    G4String limitedStr= LimitedString(fLimitedStep[num]); 
    G4cout << " " << std::setw(15) << limitedStr << " ";  
    G4cout.precision(oldPrec); 

    G4Navigator *pNav= GetNavigator( num ); 
    G4String  WorldName( "Not-Set" ); 
    if (pNav)
    {
       G4VPhysicalVolume *pWorld= pNav->GetWorldVolume(); 
       if( pWorld )
       {
           WorldName = pWorld->GetName(); 
       }
    }
    G4cout << " " << WorldName ; 
    G4cout << G4endl;
  }

  if( fVerboseLevel > 4 )
  {
    G4cout << " G4PathFinder::PrintLimited - exiting. " << G4endl;
  }
}

G4double
G4PathFinder::DoNextCurvedStep( const G4FieldTrack &initialState,
                                G4double      proposedStepLength,
                                G4VPhysicalVolume* pCurrentPhysicalVolume )
{
  const G4double toleratedRelativeError= 1.0e-10; 
  G4double minStep= kInfinity, newSafety=0.0;
  G4int numNav; 
  G4FieldTrack  fieldTrack= initialState;
  G4ThreeVector startPoint= initialState.GetPosition(); 


  G4EquationOfMotion* equationOfMotion = 
     (fpFieldPropagator->GetChordFinder()->GetIntegrationDriver()->GetStepper())
     ->GetEquationOfMotion();

  equationOfMotion->SetChargeMomentumMass( *(initialState.GetChargeState()), 
                                           initialState.GetMomentum().mag2(),
                                           initialState.GetRestMass() );
  
#ifdef G4DEBUG_PATHFINDER
  G4int prc= G4cout.precision(9);
  if( fVerboseLevel > 2 )
  {
    G4cout << " G4PathFinder::DoNextCurvedStep ****** " << G4endl;
    G4cout << " Initial value of field track is " << fieldTrack 
           << " and proposed step= " << proposedStepLength  << G4endl;
  }
#endif

  fPreStepCenterRenewed= true; // Always update PreSafety with PreStep point

  if( fNoActiveNavigators > 1 )
  { 
     // Calculate the safety values before making the step

     G4double minSafety= kInfinity, safety; 
     for( numNav=0; numNav < fNoActiveNavigators; ++numNav )
     {
        safety= fpNavigator[numNav]->ComputeSafety( startPoint, DBL_MAX, false );
        fPreSafetyValues[numNav]= safety; 
        fCurrentPreStepSafety[numNav]= safety; 
        minSafety = std::min( safety, minSafety ); 
     }

     // Save safety value, related position

     fPreSafetyLocation=  startPoint;   
     fPreSafetyMinValue=  minSafety;
     fPreStepLocation=    startPoint;
     fMinSafety_PreStepPt= minSafety;
  }

  // Allow Propagator In Field to do the hard work, calling G4MultiNavigator
  //
  minStep=  fpFieldPropagator->ComputeStep( fieldTrack,
                                            proposedStepLength,
                                            newSafety, 
                                            pCurrentPhysicalVolume );

  // fieldTrack now contains the endpoint information
  //
  fEndState= fieldTrack; 
  fMinStep=   minStep; 
  fTrueMinStep = std::min( minStep, proposedStepLength );

  if( fNoActiveNavigators== 1 )
  { 
     // Update the 'PreSafety' sphere - as any ComputeStep was called 
     // (must be done anyway in field)

     fPreSafetyValues[0]=   newSafety;
     fPreSafetyLocation= startPoint;   
     fPreSafetyMinValue= newSafety;

     // Update the current 'PreStep' point's values - mandatory
     //
     fCurrentPreStepSafety[0]= newSafety; 
     fPreStepLocation=  startPoint;
     fMinSafety_PreStepPt= newSafety;
  }

#ifdef G4DEBUG_PATHFINDER
  if( fVerboseLevel > 2 )
  {
    G4cout << "G4PathFinder::DoNextCurvedStep : " << G4endl
           << " initialState = " << initialState << G4endl
           << " and endState = " << fEndState << G4endl;
    G4cout << "G4PathFinder::DoNextCurvedStep : " 
           << " minStep = " << minStep 
           << " proposedStepLength " << proposedStepLength 
           << " safety = " << newSafety << G4endl;
  }
#endif
  G4double currentStepSize;   // = 0.0; 
  if( minStep < proposedStepLength ) // if == , then a boundary found at end ??
  {   
    // Recover the remaining information from MultiNavigator
    // especially regarding which Navigator limited the step

    G4int noLimited= 0;  //   No geometries limiting step
    for( numNav=0; numNav < fNoActiveNavigators; ++numNav )
    {
      G4double finalStep, lastPreSafety=0.0, minStepLast;
      ELimited didLimit; 
      G4bool limited; 

      finalStep=  fpMultiNavigator->ObtainFinalStep( numNav, lastPreSafety, 
                                                     minStepLast, didLimit );

      // Calculate the step for this geometry, using the 
      // final step (the only one which can differ.)

      currentStepSize = fTrueMinStep;  
      G4double diffStep= 0.0; 
      if( (minStepLast != kInfinity) )
      { 
        diffStep = (finalStep-minStepLast);
        if ( std::abs(diffStep) <= toleratedRelativeError * finalStep ) 
        {
          diffStep = 0.0;
        }
        currentStepSize += diffStep; 
      }
      fCurrentStepSize[numNav] = currentStepSize;  
      
      // TODO: could refine the way to obtain safeties for > 1 geometries
      //     - for pre step safety
      //        notify MultiNavigator about new set of sub-steps
      //        allow it to return this value in ObtainFinalStep 
      //        instead of lastPreSafety (or as well?)
      //     - for final step start (available)
      //        get final Step start from MultiNavigator
      //        and corresponding safety values
      // and/or ALSO calculate ComputeSafety at endpoint
      //     endSafety= fpNavigator[numNav]->ComputeSafety( endPoint ); 

      fLimitedStep[numNav] = didLimit; 
      fLimitTruth[numNav] = limited = (didLimit != kDoNot ); 
      if( limited ) { noLimited++; }

#ifdef G4DEBUG_PATHFINDER
      G4bool StepError= (currentStepSize < 0) 
                   || ( (minStepLast != kInfinity) && (diffStep < 0) ) ; 
      if( StepError || (fVerboseLevel > 2) )
      {
        G4String  limitedString=  LimitedString( fLimitedStep[numNav] ); 
        
        G4cout << " G4PathFinder::ComputeStep. Geometry " << numNav
               << "  step= " << fCurrentStepSize[numNav] 
               << " from final-step= " << finalStep 
               << " fTrueMinStep= " << fTrueMinStep 
               << " minStepLast= "  << minStepLast 
               << "  limited = " << (fLimitTruth[numNav] ? "YES" : " NO")
               << " ";
        G4cout << "  status = " << limitedString << " #= " << didLimit
               << G4endl;
        
        if( StepError )
        { 
          std::ostringstream message;
          message << "Incorrect calculation of step size for one navigator"
                  << G4endl
                  << "        currentStepSize = " << currentStepSize 
                  << ", diffStep= " << diffStep << G4endl
                  << "ERROR in computing step size for this navigator.";
          G4Exception("G4PathFinder::DoNextCurvedStep",
                      "GeomNav0003", FatalException, message);
        }
      }
#endif
    } // for num Navigators

    fNoGeometriesLimiting= noLimited;  // Save # processes limiting step
  } 
  else if ( (minStep == proposedStepLength)  
            || (minStep == kInfinity)  
            || ( std::abs(minStep-proposedStepLength)
               < toleratedRelativeError * proposedStepLength ) )
  { 
    // In case the step was not limited, use default responses
    //  --> all Navigators 
    // Also avoid problems in case of PathFinder using safety to optimise
    //  - it is possible that the Navigators were not called
    //    if the safety was already satisfactory.
    //    (In that case calling ObtainFinalStep gives invalid results.)

    currentStepSize= minStep;
    for( numNav=0; numNav < fNoActiveNavigators; ++numNav )
    {
      fCurrentStepSize[numNav] = minStep; 
      // Safety for endpoint ??  // Can eventuall improve it -- see TODO above
      fLimitedStep[numNav] = kDoNot; 
      fLimitTruth[numNav] = false; 
    }
    fNoGeometriesLimiting= 0;  // Save # processes limiting step
  } 
  else    //  (minStep > proposedStepLength) and not (minStep == kInfinity)
  {
    std::ostringstream message;
    message << "Incorrect calculation of step size for one navigator." << G4endl
            << "        currentStepSize = " << minStep << " is larger than "
            << " proposed StepSize = " << proposedStepLength << ".";
    G4Exception("G4PathFinder::DoNextCurvedStep()",
                "GeomNav0003", FatalException, message); 
  }

#ifdef G4DEBUG_PATHFINDER
  if( fVerboseLevel > 2 )
  {
    G4cout << " Exiting G4PathFinder::DoNextCurvedStep " << G4endl;
    PrintLimited(); 
  }
  G4cout.precision(prc); 
#endif

  return minStep; 
}


G4bool G4PathFinder::RecheckDistanceToCurrentBoundary(
                                        const G4ThreeVector &pGlobalPoint,
                                        const G4ThreeVector &pDirection,
                                        const G4double aProposedMove,
                                        G4double  *prDistance,
                                        G4double  *prNewSafety)const
{
  G4bool retval= true;
  
  if( fNoActiveNavigators > 0 )
  {
    // Calculate the safety values before making the step
    
    G4double minSafety= kInfinity;
    G4double minMove=   kInfinity;
    int numNav;
    for( numNav=0; numNav < fNoActiveNavigators; ++numNav )
    {
      G4double distance, safety;
      G4bool   moveIsOK;
      moveIsOK= fpNavigator[numNav]->RecheckDistanceToCurrentBoundary(
                                                                pGlobalPoint,
                                                                pDirection,
                                                                aProposedMove,
                                                                &distance,
                                                                &safety);
      minSafety = std::min( safety, minSafety );
      minMove   = std::min( distance, minMove );
      // The first surface encountered will determine it 
      //   - even if it is at a negative distance.
      retval &= moveIsOK;
    }
    
    *prDistance= minMove;
    if( prNewSafety ) *prNewSafety= minSafety;
  
  }else{
    retval= false;
  }

  return retval;
}



G4String& G4PathFinder::LimitedString( ELimited lim )
{
  static G4String StrDoNot("DoNot"),
                  StrUnique("Unique"),
                  StrUndefined("Undefined"),
                  StrSharedTransport("SharedTransport"),  
                  StrSharedOther("SharedOther");

  G4String* limitedStr;
  switch ( lim )
  {
     case kDoNot:  limitedStr= &StrDoNot; break;
     case kUnique: limitedStr = &StrUnique; break; 
     case kSharedTransport:  limitedStr= &StrSharedTransport; break; 
     case kSharedOther: limitedStr = &StrSharedOther; break;
     default: limitedStr = &StrUndefined; break;
  }
  return *limitedStr;
}

void G4PathFinder::PushPostSafetyToPreSafety()
{
  fPreSafetyLocation= fSafetyLocation;
  fPreSafetyMinValue= fMinSafety_atSafLocation;
  for( G4int nav=0; nav < fNoActiveNavigators; ++nav )
  {
     fPreSafetyValues[nav]= fNewSafetyComputed[nav];
  }
}
