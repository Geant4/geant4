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
// $Id: G4PathFinder.cc,v 1.7 2006-05-26 22:08:55 japost Exp $
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
       fVerboseLevel(1)
{
   fNoActiveNavigators= 0; 

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
  }
  G4cout << " PF::ComputeStep: step= " << stepNo 
	 << " nav = " << navigatorId 
	 << " trial-step-len " << proposedStepLength
	 << " from " << InitialFieldTrack.GetPosition()
	 << " dir  " << InitialFieldTrack.GetMomentumDirection()
	 << G4endl;
 
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
    // if( fVerboseLevel > 1 )
      G4cout << " current = " << currentState << G4endl;

    fCurrentStepNo = stepNo; 
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

  G4cout << " G4PathFinder returns " << fCurrentStepSize[ navigatorNo ]
	 << " for Navigator " << navigatorNo << G4endl; 

  return fCurrentStepSize[ navigatorNo ];
}

static G4TransportationManager* pTransportManager= 
       G4TransportationManager::GetTransportationManager();

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

  Locate( position, direction, false );   
  // The first location for each Navigator must be non-relative
  //   or else call ResetStackAndState() for each Navigator

  fRelocatedPoint= false; 

  if( fVerboseLevel > 1 ) 
    G4cout << " G4PathFinder::PrepareNewTrack : exiting. " << G4endl;
}

void
G4PathFinder::MovePoint()
{
  fRelocatedPoint= true;
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

  G4ThreeVector lastPosition= fEndState.GetPosition(); 
  G4ThreeVector moveVec = (position - lastPosition );
  if( (!fNewTrack) && (!fRelocatedPoint) && (moveVec.mag2() > 0.0) ){
     ReportMove( position, lastPosition, "Position" ); 
     G4Exception( "G4PathFinder::Locate", "LocateUnexpectedPoint", 
	 	  FatalException,  "Location is not where last ComputeStep ended."); 
  }

  if( fVerboseLevel > 2 ){
    G4cout << G4endl << G4endl; 
    G4cout << " G4PathFinder::Locate : entered " << G4endl;
    G4cout << " --------------------   -------" <<  G4endl;
    G4cout << "   Locating at position " << position
	   << "  with direction " << direction 
	   << "  relative= " << relative << G4endl;
  }

  for ( num=0; num< fNoActiveNavigators ; ++pNavIter,++num ) {
     //  ... who limited the step ....

     // G4Navigator pNav= *pNavIter;
     // if( fGeometryLimitedStep ) {
     if( fLimitTruth[num] ) { (*pNavIter)->SetGeometricallyLimitedStep(); }

     G4VPhysicalVolume *pLocated= 
     (*pNavIter)->LocateGlobalPointAndSetup( position, &direction,
     //*************************************//
					     relative,     // relative,
					     false);   // direction
     // Set the state related to the location
     fLocatedVolume[num] = pLocated; 
     // fStateId = kLocated; 

     // } else {
     //    (*pNavIter)->LocateGlobalPointWithinVolume( position ); 
     // }

     // Clear state related to the step
     fLimitedStep[num]   = kDoNot; 
     fCurrentStepSize[num] = 0.0;      
    
     G4cout << " Located in world " << num 
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

  if( fVerboseLevel > 2 ){
    G4cout << " G4PathFinder::Locate : exiting. " << G4endl;
    G4cout << G4endl;
  }
  fRelocatedPoint= false;
}

G4TouchableHandle 
G4PathFinder::CreateTouchableHandle( G4int navId ) const
   // Also? G4TouchableCreator& GetTouchableCreator( navId ) const; 
{
G4cout << "G4PathFinder::CreateTouchableHandle : navId = " << navId << " -- " << GetNavigator(navId) << G4endl;
  G4TouchableHistory* touchHist;
  touchHist= GetNavigator(navId) -> CreateTouchableHistory(); 
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
  G4Exception( " G4PathFinder::ComputeLinearStep is Null" );
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
G4cout << "G4PathFinder::DoNextLinearStep : Navigator [" << num << "] -- step size " << step << G4endl;
  } 
  fMinSafety= minSafety;
  fMinStep=   minStep; 

  if( fMinStep == kInfinity ){
     minStep = proposedStepLength;   //  Use this below for endpoint !!
  }
  fTrueMinStep = minStep;

  // Set the EndState
  G4ThreeVector endPosition;

G4cout << "G4PathFinder::DoNextLinearStep : initialPosition : " << initialPosition << G4endl;
  fEndState= initialState; 
  endPosition= initialPosition + minStep * initialDirection ; 
G4cout << "G4PathFinder::DoNextLinearStep : endPosition : " << endPosition << G4endl;
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
  //  static G4bool Limited[ MaxNav ]; 
  const G4int  IdTransport= 1;
  G4int noLimited=0; 
  // ELimited defaultElim= kDoNot;
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
  if( last > 0 ){ 
    if( noLimited == 1 ) {
      fLimitedStep[ last ] = kUnique; 
    }
#if 0
    else{
      // Processes limiting are not unique
      for ( num= 0; num < fNoActiveNavigators; num++ ) { 
	if( Limited[num] ) fLimitedStep[ last ] = kUnique; 
      }
    }
#endif
  }

#ifndef G4NO_VERBOSE
  this->PrintLimited();   // --> for tracing 
  if( fVerboseLevel > 2 ){
    // this->PrintLimited();  
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
  G4cout << std::setw(6) << " Step# "  << " "
	 << std::setw(6) << " NavId "  << " "
	 << std::setw(10) << " step-size " 
	 << std::setw(24) << " Limited:   bool, flag" 
         << std::setw(10) << "  World " 
	 << G4endl;  
  }
  int num;
  for ( num= 0; num < fNoActiveNavigators; num++ ) { 
    G4double stepLen = fCurrentStepSize[num]; 
    if( stepLen > fTrueMinStep ) { 
      stepLen = fTrueMinStep;     // did not limit (went as far as asked)
    }
    G4cout << std::setw(6) << fCurrentStepNo  << " " 
	   << std::setw(6) << num  << " "
	   << std::setw(10) << stepLen
	   << std::setw(8) << "         "
	   << std::setw(8) << fLimitTruth[num] << " ";
    G4String limitedStr;
    switch ( fLimitedStep[num] ) {
      case kDoNot:  limitedStr= StrDoNot; break;
      case kUnique: limitedStr = StrUnique; break; 
      case kSharedTransport:  limitedStr= StrSharedTransport; break; 
      case kSharedOther: limitedStr = StrSharedOther; break;
      default: limitedStr = StrUndefined; break;
    }
    G4cout << " " << std::setw(8) << limitedStr;  

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
