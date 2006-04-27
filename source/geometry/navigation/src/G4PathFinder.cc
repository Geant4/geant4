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
// $Id: G4PathFinder.cc,v 1.1 2006-04-27 15:28:26 japost Exp $
// GEANT4 tag $ Name:  $
// 
// class G4Navigator Implementation
//
// Original author: Paul Kent, July 95/96
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

// class G4VCurvedTrajectoryFilter;

// ********************************************************************
// Constructor
// ********************************************************************
//

static G4PathFinder*
G4PathFinder::GetInstance()
{
   static G4PathFinder* fpInstance= 0; 
   if( ! fpInstance ) {
     fpInstance= new G4PathFinder(); 
   }
   return fpInstance;
}

G4PathFinder::G4PathFinder() :
   fpActiveNavigators()
{
   // fpActiveNavigators= new std::vec<G4Navigator>;  // Null 
   fNoActiveNavigators= 0; 
   fNoTotalNavigators= 0; 

   fpTransportManager= G4TransportationManager()::GetInstance();
}

G4PathFinder::~G4PathFinder() {}

G4double 
G4PathFinder::ComputeStep( const G4FieldTrack      &InitialFieldTrack,   // Or update non-c
				 G4double           proposedStepLength,
				 G4int              navigatorId, 
				 G4int              stepNo,     // See next step / check 
				 G4double          &pNewSafety,   // for this geom 
				 G4bool            &limitedStep, 
				 G4FieldTrack      &EndState )
{
  // ---
  static G4int lastStepNo= -1;
  if( stepNo != lastStepNo ){
    // 
    G4FieldTrack currentState= InitialFieldTrack;
    DoNextCurvedStep( currentState, 
    lastStepNo= stepNo; 
  }

  // Prepare the information to return
  pNewSafety  = fNewSafety[ navId ]; 
  limitedStep = fLimitedStep[ navId ];
  EndState = fEndState; 
}

// Check and cache set of active navigators
// 
void
G4PathFinder::PrepareNewTrack( const G4ThreeVector position, const G4ThreeVector direction)
{
  // Message the G4NavigatorPanel / Dispatcher to find active navigators
  G4ThreeVector point(0.0, 0.0, 0.0); 

  std::vector<G4Navigator*>::iterator pNav= fpTransportManager-> GetActiveNavigatorsIterator();
  G4int num=0; 
  for ( pNav= fActiveNavigators::begin(); pNav != fActiveNavigators::end(); pNav++ ) {

     // Keep information in carray?
     fpNavigator[num]=  pNav;   // OR *(pNav);  ???
     // fStatus[numNav]= 
     fLimitedStep[num] = false; 
     fCurrentStepSize[num] = 0.0; 

     num++;  
  }
  fpNoActiveNavigators= num;

  Locate( position, direction ); 
}

void
G4PathFinder::Locate( const G4ThreeVector& position, const G4ThreeVector& direction)
{
  // Locate the point in each geometry
  std::vec<G4Navigator>::iterator pNav; 
  G4int num=0; 

  for ( pNav= fActiveNavigators::begin(); pNav != fActiveNavigators::end(); pNav++ ) {
     pNav->LocateGlobalPointAndSetup( position, direction, true, false); // use relative, direction
     //*****************************//
     fLimitedStep[num] = false; 
     fCurrentStepSize[num] = 0.0;      
     num++;  
  }
}

G4TouchableHistoryHandle 
G4PathFinder::CreateTouchableHandle( G4int navId ) const
   // Also? G4TouchableCreator& GetTouchableCreator( navId ) const; 
{
  return fpNavigator[navId] -> CreateTouchableHandle(); 
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

}

G4double
G4PathFinder::DoNextCurvedStep( G4FieldTrack &pFieldTrack,
				G4double      proposedStepLength
			        )
{
  // ---
  std::vec<G4Navigator>::iterator pNav; 
  G4int num=0; 
  // for ( pNav= fActiveNavigators::begin(); pNav != fActiveNavigators::end(); pNav++ ) {
  for ( num= 0; num < fNoActiveNavigators; num++ ) { 
     // Keep information in carray?
     navigator= fpNavigator[num]; 

     navigator->Compute
     // fStatus[numNav]= 
     fLimitedStep[num] = false; 
     fCurrentStepSize[num] = 0.0; 
  }
   
   
}
G4double
G4PathFinder::DoNextLinearStep( G4FieldTrack &pFieldTrack,
				G4double      proposedStepLength
			        )
{
  // ---
  std::vec<G4Navigator>::iterator pNav; 
  G4int num=0; 
  // for ( pNav= fActiveNavigators::begin(); pNav != fActiveNavigators::end(); pNav++ ) {
  for ( num= 0; num < fNoActiveNavigators; num++ ) { 
     // Keep information in carray?
     navigator= fpNavigator[num]; 

     navigator->Compute
     // fStatus[numNav]= 
     fLimitedStep[num] = false; 
     fCurrentStepSize[num] = 0.0; 
  }
   
   
}
