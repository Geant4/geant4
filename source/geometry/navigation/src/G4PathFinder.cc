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
// $Id: G4PathFinder.cc,v 1.2 2006-04-28 17:22:31 japost Exp $
// GEANT4 tag $ Name:  $
// 
// class G4PathFinder Implementation
//
// Original author:  John Apostolakis,  April 2006
//
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
  : fEndState( G4ThreeVector(), G4ThreeVector(), 0., 0., 0., 0., 0.)
{
   // fpActiveNavigators= new std::vec<G4Navigator>;  // Null 
   fNoActiveNavigators= 0; 
   // fNoNavigators= 0; 

   fMinSafety= -1.0;  // Invalid value
   fMinStep=   -1.0;  // 
   fNewTrack= false; 

   // fpNavigator= new[MaxNav] (G4Navigator*); 
}

G4PathFinder::~G4PathFinder() 
{
   // delete[] fpNavigator;
}

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
  static G4int lastStepNo= -1;
  G4int navigatorNo=-1; 
  G4cout << " G4PathFinder::ComputeStep - entered " << G4endl;
  G4cout << "   - stepNo = " << stepNo 
         << " navigatorId = " << navigatorId 
	 << " proposed step len = " << proposedStepLength
	 << G4endl;

  if( navigatorId < fNoActiveNavigators ){
    navigatorNo= navigatorId;
  } else { 
    G4cerr << " Navigator Id = " << navigatorId 
           << " No Active = " << fNoActiveNavigators << " . " << G4endl;
    G4Exception( "G4PathFinder::ComputeStep: Bad Navigator Id" ); 
  }

  if( (stepNo != lastStepNo) || fNewTrack ){
    // 
    G4FieldTrack currentState= InitialFieldTrack;
    // DoNextCurvedStep( currentState, proposedState ); 
    DoNextLinearStep( currentState, proposedStepLength ); 
    lastStepNo= stepNo; 
  }
  fNewTrack= false; 

  // Prepare the information to return
  pNewSafety  = fNewSafety[ navigatorNo ]; 
  limitedStep = fLimitedStep[ navigatorNo ];
  EndState = fEndState; 

  return fCurrentStepSize[ navigatorNo ];
}

static G4TransportationManager* pTransportManager= 
       G4TransportationManager::GetTransportationManager();

// Check and cache set of active navigators
// 
void
G4PathFinder::PrepareNewTrack( const G4ThreeVector position, 
                               const G4ThreeVector direction )
{
  G4cout << " G4PathFinder::PrepareNewTrack - entered " << G4endl;
  // static G4TransportationManager* pTransportManager= 
  //       G4TransportationManager::GetTransportationManager();

  fNewTrack= true; 

  // Message the G4NavigatorPanel / Dispatcher to find active navigators
  G4ThreeVector point(0.0, 0.0, 0.0); 

  std::vector<G4Navigator*>::iterator pNavigatorIter; 
  fNoActiveNavigators=  pTransportManager-> GetNoActiveNavigators();
  pNavigatorIter= pTransportManager-> GetActiveNavigatorsIterator();
  G4int num=0; 
  //  for ( ; pNavigatorIter<pNavigatorIter; ; ++pNavigatorIter ) { FIXME ---- TO DO 
  for( num=0; num< fNoActiveNavigators; ++pNavigatorIter,++num ) {
 
     // Keep information in carray ... for use in creating touchables at least
     fpNavigator[num] =  *pNavigatorIter;   // OR *(pNav);  ???
     // fStatus[numNav]= 
     fLimitedStep[num] = kDoNot;
     fCurrentStepSize[num] = 0.0; 
  }

  Locate( position, direction ); 

  G4cout << " G4PathFinder::PrepareNewTrack : exiting. " << G4endl;
}

void
G4PathFinder::Locate( const G4ThreeVector& position, const G4ThreeVector& direction)
{
  // Locate the point in each geometry
  std::vector<G4Navigator*>::iterator pNavIter= pTransportManager->GetActiveNavigatorsIterator(); 
  G4int num=0; 

  G4cout << " G4PathFinder::Locate : entered " << G4endl;

  // for ( pNav= fActiveNavigators::begin(); pNav != fActiveNavigators::end(); pNav++ ) {
  for ( num=0; num<fNoActiveNavigators ; ++pNavIter,++num ) {
     //  ... who limited the step ....

     (*pNavIter)->LocateGlobalPointAndSetup( position, &direction, true, false); // use relative, direction
     //*****************************//
     fLimitedStep[num] = kDoNot; 
     fCurrentStepSize[num] = 0.0;      
     num++;  
  }
  G4cout << " G4PathFinder::Locate : exiting. " << G4endl;
}

G4TouchableHandle 
G4PathFinder::CreateTouchableHandle( G4int navId ) const
   // Also? G4TouchableCreator& GetTouchableCreator( navId ) const; 
{
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
G4PathFinder::DoNextLinearStep( const G4FieldTrack &fieldTrack,
				      G4double      proposedStepLength
			      )
{
  // ---
  // G4Navigator* navigator; 
  G4double safety= 0.0, step=0.0;
  G4double minSafety= DBL_MAX, minStep= DBL_MAX;

  G4cout << " G4PathFinder::DoNextLinearStep : entered " << G4endl;
  std::vector<G4Navigator*>::iterator pNavigatorIter;

  pNavigatorIter= pTransportManager-> GetActiveNavigatorsIterator();

  G4int num=0; 
  for( num=0; num< fNoActiveNavigators; ++pNavigatorIter,++num ) {
     // navigator= this->GetNavigator(num); 
     safety= DBL_MAX;

     G4ThreeVector direction= fieldTrack.GetMomentumDirection();

     step= 
     (*pNavigatorIter)->ComputeStep( fieldTrack.GetPosition(), 
				     direction,  // &(fieldTrack.GetMomentumDirection()),
				     proposedStepLength,
				     safety ); 
     if( safety < minSafety ){ minSafety = safety; } 
     if( step < minStep ) { minStep= step; } 
     //  Later can reduce the proposed step to the latest minStep value

     fCurrentStepSize[num] = step; 
     fNewSafety[num]= safety; 
  } 
  fMinSafety= minSafety;
  fMinStep=   minStep; 

  this->WhichLimited(); 

  G4cout << " G4PathFinder::DoNextLinearStep : exits returning " << minStep << G4endl;
  return minStep;
}

void
G4PathFinder::WhichLimited()       // Flag which processes limited the step
{
  G4int num=-1, last=-1; 
  static G4bool Limited[ MaxNav ]; 
  G4int noLimited=0; 
  ELimited defaultElim= kDoNot, shared= kSharedOther; 

  G4cout << " G4PathFinder::WhichLimited - entered " << G4endl;

  // Assume that [0] is Mass / Transport
  G4bool transportLimited = (fCurrentStepSize[0] == fMinStep); 
  if( transportLimited ){ 
     shared= kSharedTransport;
  }

  for ( num= 0; num < fNoActiveNavigators; num++ ) { 
    G4bool limitedStep;

    limitedStep = ( fCurrentStepSize[num] == fMinStep ); 
    
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
    }else{
      // Processes limiting are not unique
      for ( num= 0; num < fNoActiveNavigators; num++ ) { 
	if( Limited[num] ) fLimitedStep[ last ] = kUnique; 
      }
    }
  }    
  G4cout << " G4PathFinder::WhichLimited - exiting. " << G4endl;
}
