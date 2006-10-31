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
// $Id: G4PathFinder.cc,v 1.18 2006-10-31 16:49:44 japost Exp $
// GEANT4 tag $ Name:  $
// 
// class G4PathFinder Implementation
//
// Original author:  John Apostolakis,  April 2006
// Revisions:
//  21.05.06 J.Apostolakis First implementation coworks with MassNavigator
//  23.05.06 M.Asai        Change to new navigator numbering of G4Transp..Manager
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

   pTransportManager= G4TransportationManager::GetTransportationManager();
}

G4PathFinder::~G4PathFinder() 
{
   // delete[] fpNavigator;
}

// static G4int lastStepNo= -1;

G4double 
G4PathFinder::ComputeStep( const G4FieldTrack &InitialFieldTrack, 
				 G4double     proposedStepLength,
				 G4int        navigatorId, 
				 G4int        stepNo,       // find next step 
				 G4double     &pNewSafety,  // for this geom 
				 ELimited     &limitedStep, 
				 G4FieldTrack &EndState )
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
    // 
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
      
      if( fVerboseLevel > 2 ) { 
	G4cout << " Calling PathFinder::Locate() from G4PathFinder::ComputeStep() " << G4endl;
      }

      Locate( newPosition, newDirection ); 
    }

    // DoNextCurvedStep( currentState, proposedState ); 

    DoNextLinearStep( currentState, proposedStepLength ); 
    //--------------
    fLastStepNo= stepNo; 
  }
  else{ 
    if( fVerboseLevel > 1 ){ 
      G4cout << " G4P::CS -> Not calling DoNextLinearStep: " 
	     << " stepNo= " << stepNo << " last= " << fLastStepNo 
	     << " new= " << fNewTrack << G4endl; 
    }
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

  if( fVerboseLevel > 1 ) 
    G4cout << " G4PathFinder::PrepareNewTrack - entered " << G4endl;
  // static G4TransportationManager* pTransportManager= 
  //       G4TransportationManager::GetTransportationManager();

  fNewTrack= true; 
  this->MovePoint();   // Signal further that the last status is wiped

  // Message the G4NavigatorPanel / Dispatcher to find active navigators
  G4ThreeVector point(0.0, 0.0, 0.0); 
  std::vector<G4Navigator*>::iterator pNavigatorIter; 

  fNoActiveNavigators=  pTransportManager-> GetNoActiveNavigators();
  if( fNoActiveNavigators > fMaxNav ){
    G4cerr << "Too many active Navigators (worlds). G4PathFinder fails." 
	   << G4endl;
    G4cout << " Fatal error: Transportation Manager reports "
	   << fNoActiveNavigators 
	   << " which is more than the number allowed = "
	   << fMaxNav << G4endl;
    G4Exception("G4PathFinder::PrepareNewTrack()",
		"TooManyNavigators",  FatalException,
		"Too many active Navigators / worlds"); 
  }

  pNavigatorIter= pTransportManager-> GetActiveNavigatorsIterator();
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
  std::vector<G4Navigator*>::iterator pNavIter= pTransportManager->GetActiveNavigatorsIterator(); 
  G4int num=0; 

  G4ThreeVector lastEndPosition= fEndState.GetPosition(); 
  G4ThreeVector moveVec = (position - lastEndPosition );
  G4double      moveLenSq= moveVec.mag2();
  if( (!fNewTrack) && (!fRelocatedPoint) && ( moveLenSq> 0.0) ){
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
       G4cout << " Located in world " << num << " at " << position
	      << "  used geomLimStp " << fLimitTruth[num]
	      << "  - found in volume " << pLocated ; 
       G4cout << "  name = '" ;       
       if( pLocated ){ 
	 G4cout << pLocated->GetName() << "'"; 
	 G4cout << " - CopyNo= " << pLocated->GetCopyNo(); 
       } else { 
	 G4cout <<  "Null'   Id: Not-Set "; 
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
  std::vector<G4Navigator*>::iterator pNavIter= pTransportManager->GetActiveNavigatorsIterator(); 
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
  if( (!fNewTrack) && ( longMoveEnd || longMoveSaf ) ){  
                      // For debugging, let's see when one or the other is long
    //  expected good condition is one below!!!
 // ( (!fNewTrack) && ( longMoveEnd && longMoveSaf ) ){  

     // Recompute ComputeSafety for end position
     revisedSafety= ComputeSafety(lastEndPosition); 

     G4double  distCheckRevisedEnd= 
         ( moveLenEndPosSq - revisedSafety * revisedSafety ); 
     G4bool  longMoveRevisedEnd= ( distCheckRevisedEnd > 0. ) ; 

     G4double  moveMinusSafety= 0.0; 
     G4double  moveLenEndPosition= std::sqrt( moveLenEndPosSq );
     // moveMinusSafety = moveLenEndPosition - revisedSafety; 
     moveMinusSafety = distCheckRevisedEnd
                       / (moveLenEndPosition + revisedSafety); 

#if 0
     // Temporary code to check above
     G4double moveMinus_old= ( moveLenEndPosition - revisedSafety ) ; 
     G4cout << " G4PF: New moveMinusSafety calculation " << G4endl
	    << "  new value= " << moveMinusSafety << G4endl
	    << "  old value= " << moveLenEndPosition - revisedSafety << G4endl
	    << "  difference= " << moveMinusSafety - moveMinus_old
	    << "   relative = " << (moveMinusSafety - moveMinus_old)/moveMinusSafety
	    << G4endl;
     if( fabs( moveMinusSafety - ( moveLenEndPosition - revisedSafety ) ) 
	 >  cErrorTolerance  * (moveLenEndPosition + revisedSafety) ) { 
       G4Exception( "G4PathFinder::Relocate >>> Error in move-safety calculation", 
		    "204-PathFinderMoveSafetyError", FatalException, 
		    "New method using a-b= (a^2-b^2)/(a+b) failed. ") ; 
     } 
#endif

     if( (moveMinusSafety > 0.0 ) && (revisedSafety > 0.0) ){
        // Take into account possibility of roundoff error causing
        //   this apparent move further than safety

        G4cout << " G4PF:Relocate> Ratio to revised safety is " 
	       << fabs(moveMinusSafety)/revisedSafety << G4endl;
        // 
        G4double  absMoveMinusSafety= std::fabs(moveMinusSafety);
	G4bool smallRatio= absMoveMinusSafety < kRadTolerance * revisedSafety ; 
	G4double maxCoordPos = std::max( 
				      std::max( std::fabs(position.x()), 
						std::fabs(position.y())), 
				      std::fabs(position.z()) );
	G4bool smallValue= absMoveMinusSafety < cErrorTolerance * maxCoordPos;
        if( smallRatio || smallValue ) { 
	  moveMinusSafety = 0.0; 
	  longMoveRevisedEnd = false;   // Numerical issue -- not too long!
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
	}else{
	  G4cout << " Difference of move and safety is not very small." << G4endl;
	}
     }

#if 0
     G4cout << " Needed to recalculate safety at endpoint - was potentially a problem" 
	    << G4endl;
     ReportMove( lastEndPosition, position, "Position" ); 
     G4cout << "  The endpoint finally had safety = " << revisedSafety
	    << "  vs move= "  << moveLenEndPosition << G4endl
	    << "  difference = " << moveMinusSafety << G4endl
	    << "  Squared:  safety= " << revisedSafety * revisedSafety
	    << "    move= " << moveLenEndPosSq << G4endl;
#endif

     if( longMoveSaf && longMoveRevisedEnd ){  
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
	G4cout << "       move from safety location = " << sqrt(moveLenSafSq) << G4endl;
	// moveVecSafety.mag() *** (position-fSafetyLocation).mag() << G4endl;
	G4cout << "       safety - Move-from-end= " 
	       << revisedSafety - moveLenEndPosition
	       << G4endl;
	
	
	ReportMove( lastEndPosition, position, "Position" ); 
	G4Exception( "G4PathFinder::ReLocate", "205-RelocatePointTooFar", 
	 	   FatalException,  
		  "ReLocation is further than end-safety value from step-end point (and the other-safety-value around the last-called safety 'check' point.)"); 
    }
  }

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

  for ( num=0; num< fNoActiveNavigators ; ++pNavIter,++num ) {

     //  ... none limited the step

     (*pNavIter)->LocateGlobalPointWithinVolume( position ); 
     //*************************************//

     // Clear state related to the step
     fLimitedStep[num]   = kDoNot; 
     fCurrentStepSize[num] = 0.0;      
    
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
    pNavigatorIter= pTransportManager-> GetActiveNavigatorsIterator();

    G4int num=0; 
    for( num=0; num< fNoActiveNavigators; ++pNavigatorIter,++num ) {

       G4double safety;
       safety= (*pNavigatorIter)->ComputeSafety( position ); 
  
       if( safety < minSafety ){ minSafety = safety; } 
       // fNewSafety[num]= safety; 
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
    
  return G4TouchableHandle(touchHist); 
}

  // To find the field do not forget to call
  // G4FieldManager*  FindAndSetFieldManager(G4VPhysicalVolume* pCurrentPhysVol);
  // which sets and returns the correct field manager (global or local), if any.
  // Need to call it before PropagatorInField::ComputeStep is called.
 

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

G4double
G4PathFinder::DoNextCurvedStep( const G4FieldTrack &fieldTrack,
				      G4double      proposedStepLength
		              )
{
  // ---
  G4cout << " G4PathFinder::DoNextCurvedStep" << G4endl;
  G4Exception( " G4PathFinder::DoNextCurvedStep is Unfinished" );
  G4double minStep= DBL_MAX;

  G4cout << fieldTrack << proposedStepLength  << G4endl;
  // std::vector<G4Navigator*>::iterator pNavigatorIter; 
  G4int num=0;  
  for ( num= 0; num < fNoActiveNavigators; num++ ) { 

     // navigator= this->GetNavigator(num); 

     G4double step= 0;
     // navigator->ComputeStep( ...
     // fLimitedStep[num] = kDoNot; 
     // fCurrentStepSize[num] = 0.0; 
     if( step < minStep ) { minStep= step; } 
  }   
  return minStep; 
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

  pNavigatorIter= pTransportManager-> GetActiveNavigatorsIterator();

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

  if( fVerboseLevel > 2 )
    G4cout << " G4PathFinder::WhichLimited - entered " << G4endl;

  // Assume that [0] is Mass / Transport
  G4bool transportLimited = (fCurrentStepSize[IdTransport] == fMinStep); 
  if( transportLimited ){ 
     shared= kSharedTransport;
  }

  for ( num= 0; num < fNoActiveNavigators; num++ ) { 
    G4bool limitedStep;

    G4double step= fCurrentStepSize[num]; 

    limitedStep = ( step == fMinStep ); 
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

#ifndef G4NO_VERBOSE
  if( fVerboseLevel > 1 ){
    this->PrintLimited();   // --> for tracing 
    G4cout << " G4PathFinder::WhichLimited - exiting. " << G4endl;
  }
#endif

}

void
G4PathFinder::PrintLimited()
{
  static G4String StrDoNot("DoNot"), StrUnique("Unique"), StrUndefined("Undefined"),
    StrSharedTransport("SharedTransport"),  StrSharedOther("SharedOther");
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
    G4String limitedStr;
    switch ( fLimitedStep[num] ) {
      case kDoNot:  limitedStr= StrDoNot; break;
      case kUnique: limitedStr = StrUnique; break; 
      case kSharedTransport:  limitedStr= StrSharedTransport; break; 
      case kSharedOther: limitedStr = StrSharedOther; break;
      default: limitedStr = StrUndefined; break;
    }
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

  if( fVerboseLevel > 2 )
    G4cout << " G4PathFinder::PrintLimited - exiting. " << G4endl;
}
