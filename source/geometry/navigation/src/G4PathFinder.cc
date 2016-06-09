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
// $Id: G4PathFinder.cc,v 1.25 2006/12/13 15:43:32 gunter Exp $
// GEANT4 tag $ Name:  $
// 
// class G4PathFinder Implementation
//
// Original author:  John Apostolakis,  April 2006
// Revisions:
//  21.05.06 J.Apostolakis First implementation coworks with MassNavigator
//  23.05.06 M.Asai        Change to new navigator numbering of G4Transp..Manager
//  31.05.06 J.Apostolakis New Relocate method, to relocate within volume
//  25.10.06 G.Cosmo       pTransportManager ptr a data member - singleton issue
//  31.10.06 J.Apostolakis Revised Relocate, to check proposed position vs safety
// --------------------------------------------------------------------

#include "G4PathFinder.hh"
// #include "G4ios.hh"
// #include <iomanip>

// class G4VPhysicalVolume;
// #include "G4VPhysicalVolume.hh"

class G4FieldManager;  // #include "G4FieldManager.hh"
#include "G4Navigator.hh"
#include "G4PropagatorInField.hh"
#include "G4TransportationManager.hh"
#include "G4MultiNavigator.hh"

// class G4VCurvedTrajectoryFilter;

#include <iomanip>

// ********************************************************************
// Constructor
// ********************************************************************
//

G4PathFinder*
G4PathFinder::GetInstance()
{
   static G4PathFinder* fpInstance= 0; 
   if( ! fpInstance ) {
     fpInstance= new G4PathFinder(); 
   }
   return fpInstance;
}


G4PathFinder::G4PathFinder() 
  // : fpActiveNavigators()
  : fEndState( G4ThreeVector(), G4ThreeVector(), 0., 0., 0., 0., 0.),
       fRelocatedPoint(true),
       fLastStepNo(-1), 
       fVerboseLevel(0)
{
   fpMultiNavigator= new G4MultiNavigator(); 

   fpTransportManager= G4TransportationManager::GetTransportationManager();
   fpFieldPropagator = fpTransportManager->GetPropagatorInField() ;

   fNoActiveNavigators= 0; 
   G4ThreeVector  Big3Vector( DBL_MAX, DBL_MAX, DBL_MAX ); 
   fLastLocatedPosition= Big3Vector;
   fSafetyLocation= G4ThreeVector( DBL_MAX, DBL_MAX, DBL_MAX ); 
   fPreStepLocation= Big3Vector;

   fMinSafety_PreStepPt=  -1.0; 
   fMinSafety_atSafLocation= -1.0; 
   fMinSafety= -1.0;  // Invalid value
   fMinStep=   -1.0;  // 
   fNewTrack= false; 

   G4int num;
   for( num=0; num<= fMaxNav; ++num ) {
      fpNavigator[num] =  0;   
      fLimitTruth[num] = false;
      fLimitedStep[num] = kUndefLimited;
      fCurrentStepSize[num] = -1.0; 
      fLocatedVolume[num] = 0; 
   }
   // fpNavigator= new[MaxNav] (G4Navigator*); 

}

G4PathFinder::~G4PathFinder() 
{
   // delete[] fpNavigator;
   delete fpMultiNavigator; 
}

#include "G4SafetyHelper.hh"

void
G4PathFinder::EnableParallelNavigation(G4bool enableChoice)
{
   G4Navigator *navigatorForPropagation=0, *massNavigator=0;

   massNavigator= fpTransportManager->GetNavigatorForTracking(); 
   if( enableChoice ){
      // 
      navigatorForPropagation= fpMultiNavigator; 
      G4SafetyHelper::EnableParallelNavigation(true);  // Enable Msc to use PF
   }else{
      // fpNavigator[0]; // must be mass Navigator
      navigatorForPropagation= massNavigator;       
      G4SafetyHelper::EnableParallelNavigation(false); // Disable Msc from using PF
   }
   fpFieldPropagator->SetNavigatorForPropagating(navigatorForPropagation);
}

G4double 
G4PathFinder::ComputeStep( const G4FieldTrack &InitialFieldTrack, 
				 G4double     proposedStepLength,
				 G4int        navigatorId, 
				 G4int        stepNo,       // find next step 
				 G4double     &pNewSafety,  // for this geom 
				 ELimited     &limitedStep, 
			         G4FieldTrack &EndState,
			         G4VPhysicalVolume* currentVolume)
{
  // ---
  G4int navigatorNo=-1; 

  if( fVerboseLevel > 2 ){ 
    G4cout << " -------------------------" <<  G4endl;
    G4cout << " G4PathFinder::ComputeStep - entered " << G4endl;
    G4cout << "   - stepNo = "  << std::setw(4) << stepNo  << " "
	   << " navigatorId = " << std::setw(2) << navigatorId  << " "
	   << " proposed step len = " << proposedStepLength << " "
	   << G4endl;
    G4cout << " PF::ComputeStep: step= " << stepNo 
	 << " nav = " << navigatorId 
	 << " trial-step-len " << proposedStepLength
	 << " from " << InitialFieldTrack.GetPosition()
	 << " dir  " << InitialFieldTrack.GetMomentumDirection()
	 << G4endl;
  }

  if( navigatorId <= fNoActiveNavigators ){
    navigatorNo= navigatorId;
  } else { 
    G4cerr << " Navigator Id = " << navigatorId 
           << " No Active = " << fNoActiveNavigators << " . " << G4endl;
    G4Exception( "G4PathFinder::ComputeStep: Bad Navigator Id" ); 
  }

  if( fNewTrack || (stepNo != fLastStepNo)  ){

    // This is a new track or a new step, so we must make the step
    //  ( else we can simply retrieve its results for this Navigator Id )    

    // G4cout << " initial = " << InitialFieldTrack << G4endl;
    G4FieldTrack currentState= InitialFieldTrack;
    if( fVerboseLevel > 1 )
      G4cout << " current = " << currentState << G4endl;

    fCurrentStepNo = stepNo; 

    // Check whether a process shifted the position 
    //  since the last step -- by physics processes
    G4ThreeVector newPosition = InitialFieldTrack.GetPosition();   
    G4ThreeVector moveVector= newPosition - fLastLocatedPosition; 
    G4double moveLenSq= moveVector.mag2(); 
    if( moveLenSq > kCarTolerance * kCarTolerance ){ 
       G4ThreeVector newDirection = InitialFieldTrack.GetMomentumDirection();   
       if( fVerboseLevel > 2 ) { 
          G4double moveLen= std::sqrt( moveLenSq ); 
	  G4cout << " G4PathFinder::ComputeStep : Point moved since last step " 
		 << " -- at step # = " << stepNo << G4endl
		 << " by " << moveLen 
		 << " to " << newPosition << G4endl;      
       } 

       // fRelocatedPoint= true;  //  It has moved !!
       this->MovePoint(); 

#ifdef G4VERBOSE      
       if( fVerboseLevel > 2 ) { 
	 G4cout << " Calling PathFinder::Locate() from G4PathFinder::ComputeStep() " << G4endl;
       }
#endif
       // Relocate to cope with this move -- else could abort !? 
       Locate( newPosition, newDirection ); 
    }

    // Check whether the particle have an (EM) field force exerting upon it
    //
    G4double particleCharge=  currentState.GetCharge(); 
    // G4double magDipoleMoment= currentState.GetMagneticDipoleMoment();

    G4FieldManager* fieldMgr=0;
    G4bool          fieldExertsForce = false ;
    if( (particleCharge != 0.0) ) { // || ( magDipoleMoment != 0.0 ) ) {
	fieldMgr= fpFieldPropagator->FindAndSetFieldManager( currentVolume ); 
	// Protect for case where field manager has no field (= field is zero)
	fieldExertsForce =  (fieldMgr != 0) 
	                    && (fieldMgr->GetDetectorField() != 0);
    }
    // G4cout << "G4PF::CS> Field exerts force = " << fieldExertsForce << G4endl;
    if( fieldExertsForce ) {
       DoNextCurvedStep( currentState, proposedStepLength, currentVolume ); 
       //--------------
    }else{
       DoNextLinearStep( currentState, proposedStepLength ); 
       //--------------
    }
    fLastStepNo= stepNo; 
  }
  else{
    // This is neither a new track nor a new step -- just another 
    //  client accessing information for the current track, step 
    // We will simply retrieve the results of the synchronous
    //  stepping for this Navigator Id below.
#ifdef G4VERBOSE      
    if( fVerboseLevel > 1 ){ 
      G4cout << " G4P::CS -> Not calling DoNextLinearStep: " 
	     << " stepNo= " << stepNo << " last= " << fLastStepNo 
	     << " new= " << fNewTrack << " Step already done" << G4endl; 
    }
#endif
  } 

  fNewTrack= false; 

  // Prepare the information to return
  pNewSafety  = fNewSafety[ navigatorNo ]; 
  limitedStep = fLimitedStep[ navigatorNo ];
  EndState = fEndState;  // set in DoNextLinearStep  ( right ? )
  fRelocatedPoint= false;


  if( fVerboseLevel > 0 ){ 
    G4cout << " G4PathFinder::ComputeStep returns " << fCurrentStepSize[ navigatorNo ]
	   << " for Navigator " << navigatorNo 
           << " Limited step = " << limitedStep 
           << " Safety(mm) = " << pNewSafety / mm 
	   << G4endl; 
  }

  return fCurrentStepSize[ navigatorNo ];
}

// ----------------------------------------------------------------------

void
G4PathFinder::PrepareNewTrack( const G4ThreeVector position, 
                               const G4ThreeVector direction )
{
  // Key purposes:
  //   - Check and cache set of active navigators
  //   - Reset state for new track
  G4int num=0; 

  EnableParallelNavigation(true); 
    // Switch PropagatorInField to use MultiNavigator

  if( fVerboseLevel > 1 ) 
    G4cout << " G4PathFinder::PrepareNewTrack - entered " << G4endl;
  // static G4TransportationManager* fpTransportManager= 
  //       G4TransportationManager::GetTransportationManager();

  fNewTrack= true; 
  this->MovePoint();   // Signal further that the last status is wiped

  // Message the G4NavigatorPanel / Dispatcher to find active navigators
  // G4ThreeVector point(0.0, 0.0, 0.0); 
  std::vector<G4Navigator*>::iterator pNavigatorIter; 

  fNoActiveNavigators=  fpTransportManager-> GetNoActiveNavigators();
  if( fNoActiveNavigators > fMaxNav ){
    G4cerr << "Too many active Navigators (worlds). G4PathFinder fails." << G4endl;
    G4cout << " Fatal error: Transportation Manager has "<< fNoActiveNavigators 
	   << " active navigators.  " 
	   << " This is more than the number allowed = " << fMaxNav << G4endl;
    G4Exception("G4PathFinder::PrepareNewTrack()",  "TooManyNavigators",  
		FatalException, "Too many active Navigators / worlds"); 
  }

  fpMultiNavigator->PrepareNavigators(); 
  //------------------------------------

  pNavigatorIter= fpTransportManager-> GetActiveNavigatorsIterator();
  for( num=0; num< fNoActiveNavigators; ++pNavigatorIter,++num ) {
 
     // Keep information in carray ... for creating touchables - at least
     fpNavigator[num] =  *pNavigatorIter;   
     fLimitTruth[num] = false;
     fLimitedStep[num] = kDoNot;
     fCurrentStepSize[num] = 0.0; 
     fLocatedVolume[num] = 0; 
  }

  if( fVerboseLevel > 1 ) 
  G4cout << " Calling PathFinder::Locate() from G4PathFinder::PrepareNewTrack() " << G4endl;

  Locate( position, direction, false );   
  // The first location for each Navigator must be non-relative
  //   or else call ResetStackAndState() for each Navigator

  fRelocatedPoint= false; 

  if( fVerboseLevel > 0 ) {
    G4cout << " G4PathFinder::PrepareNewTrack : exiting. " << G4endl;
  }
}

static 
void ReportMove( G4ThreeVector OldVector, G4ThreeVector NewVector, G4String Quantity )
{
    G4ThreeVector moveVec = ( NewVector - OldVector );

    G4cerr << G4endl
	   << "**************************************************************" << G4endl;
    G4cerr << "Endpoint has moved between value returned by ComputeStep"
	   << " and call to Locate. " << G4endl
	   << "Change of " << Quantity << " is " << moveVec.mag() / mm << " mm long, "
	   << " and its vector is " << (1.0/mm) * moveVec << " mm " << G4endl
	   << "Endpoint of ComputeStep was " << OldVector
	   << " and current position to locate is " << NewVector << G4endl;
}

void
G4PathFinder::Locate( const   G4ThreeVector& position, 
		      const   G4ThreeVector& direction,
		      G4bool  relative)
{
  // Locate the point in each geometry
  std::vector<G4Navigator*>::iterator pNavIter= fpTransportManager->GetActiveNavigatorsIterator(); 
  G4int num=0; 

  G4ThreeVector lastEndPosition= fEndState.GetPosition(); 
  G4ThreeVector moveVec = (position - lastEndPosition );
  G4double      moveLenSq= moveVec.mag2();
  if( (!fNewTrack) && (!fRelocatedPoint) && ( moveLenSq> kCarTolerance*kCarTolerance ) ){   // ( moveLenSq> 0.0) ){
     ReportMove( position, lastEndPosition, "Position" ); 
     G4Exception( "G4PathFinder::Locate", "201-LocateUnexpectedPoint", 
	 	  JustWarning,   // FatalException,  
		  "Location is not where last ComputeStep ended."); 
  }
  fLastLocatedPosition= position; 

  if( fVerboseLevel > 2 ){
    G4cout << G4endl; 
    G4cout << " G4PathFinder::Locate : entered " << G4endl;
    G4cout << " --------------------   -------" <<  G4endl;
    G4cout << "   Locating at position " << position << "  with direction " << direction 
	   << "  relative= " << relative << G4endl;
    if ( (fVerboseLevel > 1) || ( moveLenSq > 0.0) ){ 
       G4cout << "  lastEndPosition = " << lastEndPosition
	      << "  moveVec = " << moveVec
	      << "  newTr = " << fNewTrack 
	      << "  relocated = " << fRelocatedPoint << G4endl;
    }
  }

  if( fVerboseLevel > 2 ) { 
    G4cout << " Located at " << position ; 
    if( fNoActiveNavigators > 1 )  G4cout << G4endl;
  }

  for ( num=0; num< fNoActiveNavigators ; ++pNavIter,++num ) {
     //  ... who limited the step ....

     // G4Navigator pNav= *pNavIter;
     // if( fGeometryLimitedStep ) {
     if( fLimitTruth[num] ) { (*pNavIter)->SetGeometricallyLimitedStep(); }

     G4VPhysicalVolume *pLocated= 
     (*pNavIter)->LocateGlobalPointAndSetup( position, &direction,
     //*************************************//
					     relative,  
					     false);   
     // Set the state related to the location
     fLocatedVolume[num] = pLocated; 

     // Clear state related to the step
     fLimitedStep[num]   = kDoNot; 
     fCurrentStepSize[num] = 0.0;      
    
     if( fVerboseLevel > 2 ){
       G4cout << " - In world " << num 
	      << "  geomLimStp= " << fLimitTruth[num]
	      << "  gives volume= " << pLocated ; 
       if( pLocated ){ 
	 G4cout << "  name = '" << pLocated->GetName() << "'"; 
	 G4cout << " - CopyNo= " << pLocated->GetCopyNo(); 
       } 
       G4cout  << G4endl; 
     }
  } // ending for (num= ....

  if( fVerboseLevel > 2 ){
    G4cout << " G4PathFinder::Locate : exiting. " << G4endl << G4endl; 
  }
  fRelocatedPoint= false;
}

void
G4PathFinder::ReLocate( const   G4ThreeVector& position )
  //	const   G4ThreeVector& direction, G4bool relative  )
{
  // Locate the point in each geometry
  std::vector<G4Navigator*>::iterator pNavIter= fpTransportManager->GetActiveNavigatorsIterator(); 
  const G4double cErrorTolerance=1e-12;   
    // Maximum relative error from roundoff of arithmetic 
  G4int num=0; 

  // Check that this relocation does not violate safety
  //   - at endpoint (computed from start point) AND
  //   - at last safety location  (likely just called)

  G4ThreeVector lastEndPosition= fEndState.GetPosition(); 
  // calculate end-point safety 
  G4double      DistanceStartEnd= (lastEndPosition - fPreStepLocation).mag();
  G4double      endPointSafety_raw = fMinSafety_PreStepPt 
                                    - DistanceStartEnd; 
                                 // - fTrueMinStep; // OK for linear only 
  G4double      endPointSafety_Est1 = std::max( 0.0,  endPointSafety_raw ); 

  //   and check move from endpoint against this endpoint safety
  G4ThreeVector moveVecEndPos  = (position - lastEndPosition );
  G4double      moveLenEndPosSq = moveVecEndPos.mag2(); 

  // Check that move from endpoint of last step is within safety
  //    -- or check against last location or relocation ?? 
  // G4ThreeVector moveVecStepStart = (position - fLastPrestepPosition );
  // G4double      moveLenStartSq= moveVecStepStart.mag2(); 
  G4ThreeVector moveVecSafety=     (position - fSafetyLocation ); 
  G4double      moveLenSafSq=   moveVecSafety.mag2();

  // G4double distCheckStart_sq= ( moveLenStartSq - fMinSafety*fMinSafety ); 
  G4double distCheckEnd_sq= ( moveLenEndPosSq - endPointSafety_Est1 
			                       *endPointSafety_Est1 ); 
  G4double distCheckSaf_sq=   ( moveLenSafSq -  fMinSafety_atSafLocation
			                       *fMinSafety_atSafLocation ); 
  // G4bool longMoveStart = distCheckStart_sq > 0.0; 
  G4bool longMoveEnd = distCheckEnd_sq > 0.0; 
  G4bool longMoveSaf = distCheckSaf_sq > 0.0; 

  G4double revisedSafety= 0.0; 
  if( (!fNewTrack) && ( longMoveEnd && longMoveSaf ) ){  
       // Used to use ( longMoveEnd || longMoveSaf ) for extra checking

     // Recompute ComputeSafety for end position
     revisedSafety= ComputeSafety(lastEndPosition); 

     G4double  distCheckRevisedEnd= 
         ( moveLenEndPosSq - revisedSafety * revisedSafety ); 
     G4bool  longMoveRevisedEnd=  ( distCheckRevisedEnd > 0. ) ; 

     G4double  moveMinusSafety= 0.0; 
     G4double  moveLenEndPosition= std::sqrt( moveLenEndPosSq );
     moveMinusSafety = moveLenEndPosition - revisedSafety; 
     // moveMinusSafety = distCheckRevisedEnd / (moveLenEndPosition + revisedSafety); 

     if( longMoveRevisedEnd && (moveMinusSafety > 0.0 ) && (revisedSafety > 0.0) ){
        // Take into account possibility of roundoff error causing
        //   this apparent move further than safety

        if( fVerboseLevel > 0 ) 
	   G4cout << " G4PF:Relocate> Ratio to revised safety is " 
		  << std::fabs(moveMinusSafety)/revisedSafety << G4endl;
        // 
        G4double  absMoveMinusSafety= std::fabs(moveMinusSafety);
	G4bool smallRatio= absMoveMinusSafety < kRadTolerance * revisedSafety ; 
	G4double maxCoordPos = std::max( 
				      std::max( std::fabs(position.x()), 
						std::fabs(position.y())), 
				      std::fabs(position.z()) );
	G4bool smallValue= absMoveMinusSafety < cErrorTolerance * maxCoordPos;
        if( ! (smallRatio || smallValue) ) { 
	  G4cout << " G4PF:Relocate> Ratio to revised safety is " 
 	         << std::fabs(moveMinusSafety)/revisedSafety << G4endl;
	  G4cout << " Difference of move and safety is not very small." << G4endl;
	}else{
	  moveMinusSafety = 0.0; 
	  longMoveRevisedEnd = false;   // Numerical issue -- not too long!
#ifdef G4DEBUG_PATHFINDER
	  G4cout << " Difference of move and safety is very small in magnitude, " 
		 << absMoveMinusSafety << G4endl;
	  if( smallRatio ) {
	    G4cout << " ratio to safety " << revisedSafety 
		   << " is " <<  absMoveMinusSafety / revisedSafety
		   << "smaller than " << kRadTolerance << " of safety ";
	  }else{
	    G4cout << " as fraction " << absMoveMinusSafety / maxCoordPos 
		   << " of position vector max-coord " << maxCoordPos
		   << " smaller than " << cErrorTolerance ;
	  }
          G4cout << " -- reset moveMinusSafety to " << moveMinusSafety << G4endl;
#endif
	}
     }

     if ( longMoveEnd && longMoveSaf && longMoveRevisedEnd && (moveMinusSafety>0.0)) { 
        // if( (moveMinusSafety>0.0) ){   // Eventually ? 
        G4int oldPrec= G4cout.precision(9); 
        G4cout << " Problem in G4PathFinder::Relocate() " << G4endl;
        G4cout << " Moved from last endpoint by " << moveLenEndPosition 
	  // << " Moved from last located by " << std::sqrt(moveLenStartSq) 
	       << " compared to end safety (from preStep point) = " 
	       << endPointSafety_Est1 << G4endl; 

	G4cout << "  --> last PreStep Location was " << fPreStepLocation << G4endl;
	G4cout << "       safety value =  " << fMinSafety_PreStepPt << G4endl;
	G4cout << "  --> last EndStep Location was " << lastEndPosition << G4endl;
	G4cout << "       safety value =  " << endPointSafety_Est1 
	       << " raw-value = " << endPointSafety_raw << G4endl;
	G4cout << "  --> Calling again at this endpoint, we get "
	       <<  revisedSafety << " as safety value."  << G4endl;
	G4cout << "  --> last position for safety " << fSafetyLocation << G4endl;
	G4cout << "       its safety value =  " << fMinSafety_atSafLocation << G4endl;
	G4cout << "       move from safety location = " << std::sqrt(moveLenSafSq) << G4endl;
	// moveVecSafety.mag() *** (position-fSafetyLocation).mag() << G4endl;
	G4cout << "       safety - Move-from-end= " 
	       << revisedSafety - moveLenEndPosition << " (negative is Bad.)" << G4endl;

	G4cout << " Debub:  distCheckRevisedEnd = " << distCheckRevisedEnd << G4endl;

	ReportMove( lastEndPosition, position, "Position" ); 
	G4Exception( "G4PathFinder::ReLocate", "205-RelocatePointTooFar", 
	 	   FatalException,  
		  "ReLocation is further than end-safety value from step-end point (and the other-safety-value around the last-called safety 'check' point.)"); 
	G4cout.precision(oldPrec); 
    }
  }

#ifdef G4VERBOSE
  if( fVerboseLevel > 2 ){
    G4cout << G4endl; 
    G4cout << " G4PathFinder::ReLocate : entered " << G4endl;
    G4cout << " ----------------------   -------" <<  G4endl;
    G4cout << "  *Re*Locating at position " << position  << G4endl; 
      // << "  with direction " << direction 
      // << "  relative= " << relative << G4endl;
    if ( (fVerboseLevel > -1) || ( moveLenEndPosSq > 0.0) ){ 
       G4cout << "  lastEndPosition = " << lastEndPosition
	      << "  moveVec from step-end = " << moveVecEndPos
	      << "  is new Track = " << fNewTrack 
	      << "  relocated = " << fRelocatedPoint << G4endl;
    }
  }
#endif

  for ( num=0; num< fNoActiveNavigators ; ++pNavIter,++num ) {

     //  ... none limited the step

     (*pNavIter)->LocateGlobalPointWithinVolume( position ); 
     //*************************************//

     // Clear state related to the step
     fLimitedStep[num]   = kDoNot; 
     fCurrentStepSize[num] = 0.0;      
     fLimitTruth[num] = false;   
     // G4cout << " ReLocated in world " << num << " at " << position << G4endl;
  }

  fLastLocatedPosition= position; 
  fRelocatedPoint= false;

  if( fVerboseLevel > 2 ){
    G4cout << " G4PathFinder::ReLocate : exiting " 
	   << "  at position " << fLastLocatedPosition << G4endl;
    G4cout << G4endl;
  }

}

// -----------------------------------------------------------------------------

G4double  G4PathFinder::ComputeSafety( const G4ThreeVector& position )
     // Recompute safety for the relevant point
{
    G4double minSafety= DBL_MAX; 
    // G4cout << " G4PathFinder::ComputeSafety - called at " << position << G4endl;
  
    std::vector<G4Navigator*>::iterator pNavigatorIter;
    pNavigatorIter= fpTransportManager-> GetActiveNavigatorsIterator();

    G4int num=0; 
    for( num=0; num< fNoActiveNavigators; ++pNavigatorIter,++num ) {

       G4double safety;
       safety= (*pNavigatorIter)->ComputeSafety( position ); 
  
       if( safety < minSafety ){ minSafety = safety; } 
       // fNewSafety[num]= safety; 
  
       // G4cout << " Navigator # " << num << " gives safety = " << safety << G4endl;
  } 

    fSafetyLocation= position;
    fMinSafety_atSafLocation = minSafety;

    if( fVerboseLevel > 1 ) { 
      G4cout << " G4PathFinder::ComputeSafety - returns " 
	     << minSafety << " at location " << position 
	     << G4endl;
    }
    return minSafety; 
}


// -----------------------------------------------------------------------------

G4TouchableHandle 
G4PathFinder::CreateTouchableHandle( G4int navId ) const
   // Also? G4TouchableCreator& GetTouchableCreator( navId ) const; 
{
  if( fVerboseLevel > 2 ){
    G4cout << "G4PathFinder::CreateTouchableHandle : navId = " << navId << " -- " << GetNavigator(navId) << G4endl;
  }

  G4TouchableHistory* touchHist;
  touchHist= GetNavigator(navId) -> CreateTouchableHistory(); 

  // return G4TouchableHandle(touchHist); 
     //-------->  Problem with out of world !!

  // G4TouchableHistory* touchHist= new G4TouchableHistory(); 

  G4VPhysicalVolume* locatedVolume= fLocatedVolume[navId]; 
  if( locatedVolume == 0 )
     {
       // Workaround to ensure that the touchable is fixed !! // TODO: fix
       touchHist->UpdateYourself( locatedVolume, 
				  touchHist->GetHistory() );
     }
 
  if( fVerboseLevel > 2 ) {   
    G4String VolumeName("None"); 
    if( locatedVolume ) { VolumeName= locatedVolume->GetName(); }
    G4cout << " Touchable History created at address " << touchHist
	   << "  volume = " << locatedVolume << " name= " << VolumeName
	   << G4endl;
  }

  return G4TouchableHandle(touchHist); 
}

G4double
G4PathFinder::DoNextLinearStep( const G4FieldTrack &initialState,
				      G4double      proposedStepLength
			      )
{
  // ---
  // G4Navigator* navigator; 
  G4double safety= 0.0, step=0.0;
  G4double minSafety= DBL_MAX, minStep= DBL_MAX;

  if( fVerboseLevel > 2 ){
    G4cout << " G4PathFinder::DoNextLinearStep : entered " << G4endl;
    G4cout << "   Input field track= " << initialState << G4endl;
    G4cout << "   Requested step= " << proposedStepLength << G4endl;
  }

  std::vector<G4Navigator*>::iterator pNavigatorIter;

  pNavigatorIter= fpTransportManager-> GetActiveNavigatorsIterator();

  G4ThreeVector initialPosition= initialState.GetPosition(); 
  G4ThreeVector initialDirection= initialState.GetMomentumDirection();

  G4int num=0; 
  for( num=0; num< fNoActiveNavigators; ++pNavigatorIter,++num ) {
     // navigator= this->GetNavigator(num); 
     safety= DBL_MAX;

     step= 
     (*pNavigatorIter)->ComputeStep( initialPosition, 
				     initialDirection,
				     proposedStepLength,
				     safety ); 
     if( safety < minSafety ){ minSafety = safety; } 
     if( step < minStep ) { minStep= step; } 
     //  Later can reduce the proposed step to the latest minStep value

     // if( step == kInfinity ) { step = proposedStepLength; }
     fCurrentStepSize[num] = step; 
     fNewSafety[num]= safety; 

     if( fVerboseLevel > 2 ){
       G4cout << "G4PathFinder::DoNextLinearStep : Navigator [" << num << "] -- step size " << step << G4endl;
     }
  } 

  // Save safety value, related position
  fPreStepLocation=     initialPosition; 
  fMinSafety_PreStepPt= minSafety;

  // Also store in simple 'safety' status ?
  // fMinSafety= minSafety;
  // fMinSafety_atSafLocation = minSafety;

  fMinStep=   minStep; 

  if( fMinStep == kInfinity ){
     minStep = proposedStepLength;   //  Use this below for endpoint !!
  }
  fTrueMinStep = minStep;

  // Set the EndState
  G4ThreeVector endPosition;

  fEndState= initialState; 
  endPosition= initialPosition + minStep * initialDirection ; 

  if( fVerboseLevel > 1 ){
    G4cout << "G4PathFinder::DoNextLinearStep : "
	   << " initialPosition = " << initialPosition 
	   << " and endPosition = " << endPosition<< G4endl;
  }

  fEndState.SetPosition( endPosition ); 
  fEndState.SetProperTimeOfFlight( -1.000 );   // Not defined YET
  // fEndState.SetMomentum( initialState.GetMomentum ); 
  this->WhichLimited(); 

  if( fVerboseLevel > 2 ){
    G4cout << " G4PathFinder::DoNextLinearStep : exits returning " << minStep << G4endl;
    G4cout << "   Endpoint values = " << fEndState << G4endl;
    G4cout << G4endl;
  }

  return minStep;
}

void
G4PathFinder::WhichLimited()       // Flag which processes limited the step
{
  G4int num=-1, last=-1; 
  const G4int  IdTransport= 0;  // Id of Mass Navigator !!
  G4int noLimited=0; 
  ELimited shared= kSharedOther; 

  if( fVerboseLevel > 4 )
    G4cout << " G4PathFinder::WhichLimited - entered " << G4endl;

  // Assume that [IdTransport] is Mass / Transport
  G4bool transportLimited = (fCurrentStepSize[IdTransport] == fMinStep)
                           && ( fMinStep!= kInfinity) ; 
  if( transportLimited ){ 
     shared= kSharedTransport;
  }

  for ( num= 0; num < fNoActiveNavigators; num++ ) { 
    G4bool limitedStep;

    G4double step= fCurrentStepSize[num]; 

    // limitedStep = ( step == fMinStep ); 
    limitedStep = ( step == fMinStep ) && ( step != kInfinity); 
    // if( step == kInfinity ) { fCurrentStepSize[num] = proposedStepLength; }
   
    fLimitTruth[ num ] = limitedStep; 
    if( limitedStep ) {
      noLimited++;  
      fLimitedStep[num] = shared;
      last= num; 
    }else{
      fLimitedStep[num] = kDoNot;
    }
  }
  if( (last > -1) && (noLimited == 1 ) ){
    fLimitedStep[ last ] = kUnique; 
  }

#ifdef G4VERBOSE
  if( fVerboseLevel > 1 ){
    this->PrintLimited();   // --> for tracing 
    if( fVerboseLevel > 4 )
      G4cout << " G4PathFinder::WhichLimited - exiting. " << G4endl;
  }
#endif

}

void
G4PathFinder::PrintLimited()
{
  G4String& LimitedString( ELimited lim ); 

  // Report results -- for checking   
  G4cout << "G4PathFinder::PrintLimited reports: " ; 
  G4cout << "  Minimum step (true)= " << fTrueMinStep 
	 << "  reported min = " << fMinStep 
	 << G4endl; 
  if(  (fCurrentStepNo <= 2) || (fVerboseLevel>=2) ) {
  G4cout << std::setw(5) << " Step#"  << " "
	 << std::setw(5) << " NavId"  << " "
	 << std::setw(12) << " step-size " << " "
	 << std::setw(12) << " raw-size "  << " "
	 << std::setw(12) << " pre-safety " << " " 
	 << std::setw(15) << " Limited / flag"  << " "
         << std::setw(15) << "  World "  << " "
	 << G4endl;  
  }
  int num;
  for ( num= 0; num < fNoActiveNavigators; num++ ) { 
    G4double rawStep = fCurrentStepSize[num]; 
    G4double stepLen = fCurrentStepSize[num]; 
    if( stepLen > fTrueMinStep ) { 
      stepLen = fTrueMinStep;     // did not limit (went as far as asked)
    }
    G4int oldPrec= G4cout.precision(9); 
    // const char *BooleanValue[2] = { " NO", "YES" } ; 
    G4cout << std::setw(5) << fCurrentStepNo  << " " 
	   << std::setw(5) << num  << " "
	   << std::setw(12) << stepLen << " "
	   << std::setw(12) << rawStep << " "
	   << std::setw(12) << fNewSafety[num] << " "
	   << std::setw(5) << (fLimitTruth[num] ? "YES" : " NO") << " ";
    G4String limitedStr= LimitedString(fLimitedStep[num]); 
    G4cout << " " << std::setw(15) << limitedStr << " ";  
    G4cout.precision(oldPrec); 

    G4Navigator *pNav= GetNavigator( num ); 
    G4String  WorldName( "Not-Set" ); 
    if (pNav) {
       G4VPhysicalVolume *pWorld= pNav->GetWorldVolume(); 
       if( pWorld ) {
 	  WorldName = pWorld->GetName(); 
       }
    }
    G4cout << " " << WorldName ; 
    G4cout << G4endl;
  }

  if( fVerboseLevel > 4 )
    G4cout << " G4PathFinder::PrintLimited - exiting. " << G4endl;
}


G4double
G4PathFinder::DoNextCurvedStep( const G4FieldTrack &initialState,
				G4double      proposedStepLength,
				G4VPhysicalVolume* pCurrentPhysicalVolume )
{
  const G4double toleratedRelativeError= 1.0e-10; 
  if( fVerboseLevel > 2 )
    G4cout << " G4PathFinder::DoNextCurvedStep ****** " << G4endl;
  G4double minStep= DBL_MAX, newSafety=0.0;
  G4int numNav; 
  G4FieldTrack  fieldTrack= initialState;

  if( fVerboseLevel > 2 )
    G4cout << " Initial value of field track is " << fieldTrack 
	   << " and proposed step= " << proposedStepLength  << G4endl;

  // Allow Propagator In Field to do the hard work, calling G4MultiNavigator
  minStep=  fpFieldPropagator->ComputeStep( fieldTrack,
					    proposedStepLength,
					    newSafety, 
					    pCurrentPhysicalVolume );
  // fieldTrack now contains the endpoint information
  fEndState= fieldTrack; 
  fMinStep=   minStep; 
  fTrueMinStep = std::min( minStep, proposedStepLength );

  if( fVerboseLevel > 2 )
    G4cout << "G4PathFinder::DoNextCurvedStep : " << G4endl
	   << " initialState = " << initialState << G4endl
	   << " and endState = " << fEndState << G4endl;

  G4double currentStepSize = 0;  
  G4cout << "G4PathFinder::DoNextCurvedStep : " 
	 << " minStep = " << minStep 
	 << " proposedStepLength " << proposedStepLength 
	 << " safety = " << newSafety << G4endl;
  if( minStep < proposedStepLength ) {   // if == , then a boundary found at end ??

    // Recover the remaining information from MultiNavigator
    //   especially regarding which Navigator limited the step
    for( numNav=0; numNav < fNoActiveNavigators; ++numNav ) {
      G4double finalStep, lastPreSafety=0.0, minStepLast;
      ELimited didLimit; 

      finalStep=  fpMultiNavigator->ObtainFinalStep( numNav, lastPreSafety, 
						     minStepLast, didLimit );
      currentStepSize = fTrueMinStep;  
      G4double diffStep= 0.0; 
      if( (minStepLast != kInfinity) ){ 
	diffStep = (finalStep-minStepLast);
	if ( std::abs(diffStep) <= toleratedRelativeError * finalStep ) 
	  { diffStep = 0.0; } 
	currentStepSize += diffStep; 
      }
      fCurrentStepSize[numNav] = currentStepSize;  
      
      // if( (currentStepSize < 0) || (diffStep < 0) ) {  }
      
      // fNewSafety[numNav]= preSafety;
      fNewSafety[numNav]= 0.0;  
    // TODO: improve this safety !!
      // Currently would need to call ComputeSafety at start or endpoint
      //     endSafety= fpNavigator[numNav]->ComputeSafety( endPoint ); 
      //     ...
      //   else 
      //     - for pre step safety
      //        notify MultiNavigator about new set of sub-steps
      //        allow it to return this value in ObtainFinalStep 
      //        instead of lastPreSafety (or as well?)
      //     - for final step start (available)
      //        get final Step start from MultiNavigator
      fLimitedStep[numNav] = didLimit; 
      fLimitTruth[numNav] = (didLimit != kDoNot ); 
      
      G4bool StepError= (currentStepSize < 0) 
                   || ( (minStepLast != kInfinity) && (diffStep < 0) ) ; 
      if( StepError || (fVerboseLevel > 2) ){
	
	G4String& LimitedString( ELimited lim ); 
	G4String  limitedString=  LimitedString( fLimitedStep[numNav] ); 
	
	G4cout << " G4PathFinder::ComputeStep. Geometry " << numNav << "  step= " << fCurrentStepSize[numNav] 
	       << " from final-step= " << finalStep 
	       << " fTrueMinStep= " << fTrueMinStep 
	       << " minStepLast= "  << minStepLast 
	       << "  limited = " << (fLimitTruth[numNav] ? "YES" : " NO") << " ";
	G4cout << "  status = " << limitedString << " #= " << didLimit << G4endl;
	
	if( StepError ) { 
          G4cerr << " currentStepSize = " << currentStepSize 
                 << " diffStep= " << diffStep << G4endl;
 	  G4cerr << "ERROR in computing step size for this navigator." << G4endl;
	  G4Exception( "G4PathFinder::DoNextCurvedStep", "207-StepGoingBackwards", 
		       FatalException,  
		       "Incorrect calculation of step size for one navigator" ); 
	} 
      }

    } // for num Navigators

  } 
  // else if ( minStep > proposedStepLength ) { // ie minStep == kInfinity
  else if ( (minStep == proposedStepLength)  
	    || (minStep == kInfinity)  
	    || ( std::abs(minStep-proposedStepLength) < toleratedRelativeError * proposedStepLength )
	  ) { 
    // In case the step was not limited, use default responses
    //  --> all Navigators 
    // Also avoid problems in case of PathFinder using safety to optimise
    //  - it is possible that the Navigators were not called
    //    if the safety was already satisfactory.
    //    (In that case calling ObtainFinalStep gives invalid results.)
    currentStepSize= minStep;  
    for( numNav=0; numNav < fNoActiveNavigators; ++numNav ) {
      fCurrentStepSize[numNav] = minStep; 
      fNewSafety[numNav]= 0.0;  
      // Improve it -- see TODO above
      fLimitedStep[numNav] = kDoNot; 
      fLimitTruth[numNav] = false; 
    }
  } 
  else{   //  (minStep > proposedStepLength) and not (minStep == kInfinity)
    G4cerr << G4endl;
    G4cerr << "ERROR in G4PathFinder::DoNextCurvedStep "
 	   << " currentStepSize = " << minStep << " is larger than "
	   << " proposed StepSize = " << proposedStepLength 
	   << G4endl;
    G4Exception( "G4PathFinder::DoNextCurvedStep", "208-StepLongerThanRequested", 
		 FatalException,  
		 "Incorrect calculation of step size for one navigator" ); 
  }

  if( fVerboseLevel > 2 ){
    G4cout << " Exiting G4PathFinder::DoNextCurvedStep " << G4endl;
    PrintLimited(); 
  }

  return minStep; 
}

static G4String StrDoNot("DoNot"), StrUnique("Unique"), 
  StrUndefined("Undefined"), StrSharedTransport("SharedTransport"),  
  StrSharedOther("SharedOther");

G4String& LimitedString( ELimited lim )
{
  G4String* limitedStr;
  switch ( lim ) {
     case kDoNot:  limitedStr= &StrDoNot; break;
     case kUnique: limitedStr = &StrUnique; break; 
     case kSharedTransport:  limitedStr= &StrSharedTransport; break; 
     case kSharedOther: limitedStr = &StrSharedOther; break;
     default: limitedStr = &StrUndefined; break;
  }

  return *limitedStr;
}


#if 0 
// Potential extension ..... ?? 
     // Return relevant step
     //  When no field exists or the particle has no charge or EM moment

G4double 
G4PathFinder::ComputeLinearStep(const G4ThreeVector &pGlobalPoint,
                              const G4ThreeVector &pDirection,
                              G4double pCurrentProposedStepLength,
                              G4double  &pNewSafety,
                              G4bool    &limitedStep, 
                              G4int     stepNo,       // See next step / check 
                              G4int     navId ) 
{
  G4cout << pGlobalPoint << pDirection << pCurrentProposedStepLength
	 << pNewSafety << limitedStep << stepNo << navId << G4endl; 

  G4cout << " G4PathFinder::ComputeLinearStep" << G4endl;
  G4Exception( " G4PathFinder::ComputeLinearStep is Null",  "203-No Method",
	       FatalException,  "G4PathFinder::ComputeLinearStep is Null" );
  return 0.500 * pCurrentProposedStepLength;
}
#endif
