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
// $Id: G4MultiNavigator.cc,v 1.4 2006/11/14 15:41:56 japost Exp $
// GEANT4 tag $ Name:  $
// 
// class G4PathFinder Implementation
//
// Author:  John Apostolakis, November 2006
// --------------------------------------------------------------------

#include "G4MultiNavigator.hh"

class G4FieldManager;

#include "G4Navigator.hh"
#include "G4PropagatorInField.hh"
#include "G4TransportationManager.hh"

#include <iomanip>

// ********************************************************************
// Constructor
// ********************************************************************
//
G4MultiNavigator::G4MultiNavigator() 
  // : fpActiveNavigators()
  : G4Navigator(), 
    fVerboseLevel(1)
{
   fNoActiveNavigators= 0; 
   G4ThreeVector  Big3Vector( DBL_MAX, DBL_MAX, DBL_MAX ); 
   fLastLocatedPosition= Big3Vector;
   fSafetyLocation=  Big3Vector;
   fPreStepLocation= Big3Vector;

   fMinSafety_PreStepPt=  -1.0; 
   fMinSafety_atSafLocation= -1.0; 
   fMinSafety= -DBL_MAX;  
   fMinStep=   -DBL_MAX;  
   // fNewTrack= false; 

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

   // EndState = G4FieldTrack( G4ThreeVector(), G4ThreeVector(), 0., 0., 0., 0., 0.) ); 
   // fRelocatedPoint(
   // fLastStepNo= -1;  

   G4Navigator* massNav= pTransportManager->GetNavigatorForTracking();
   if( massNav ) { 
     G4VPhysicalVolume* pWorld= massNav->GetWorldVolume(); 
     if( pWorld ) { 
       this->SetWorldVolume( pWorld ); 
       fLastMassWorld= pWorld; 
     }
   }
}

G4MultiNavigator::~G4MultiNavigator() 
{
   // delete[] fpNavigator;
}

// static G4int lastStepNo= -1;
  // To find the field do not forget to call
  // G4FieldManager*  FindAndSetFieldManager(G4VPhysicalVolume* pCurrentPhysVol);
  // which sets and returns the correct field manager (global or local), if any.
  // Need to call it before PropagatorInField::ComputeStep is called.

G4double G4MultiNavigator::ComputeStep(const G4ThreeVector &pGlobalPoint,
                     const G4ThreeVector &pDirection,
                     const G4double      proposedStepLength,
                           G4double      &pNewSafety)
{
  G4double safety= 0.0, step=0.0;
  G4double minSafety= DBL_MAX, minStep= DBL_MAX;

  if( fVerboseLevel > 2 ){
    G4cout << " G4MultiNavigator::ComputeStep : entered " << G4endl;
    G4cout << "   Input position= " << pGlobalPoint
           << "   direction= "      << pDirection         << G4endl;
    G4cout << "   Requested step= " << proposedStepLength << G4endl;
  }

  std::vector<G4Navigator*>::iterator pNavigatorIter;

  pNavigatorIter= pTransportManager-> GetActiveNavigatorsIterator();

  G4ThreeVector initialPosition=  pGlobalPoint;
  G4ThreeVector initialDirection= pDirection;

  G4int num=0; 
  for( num=0; num< fNoActiveNavigators; ++pNavigatorIter,++num ) {
     safety= DBL_MAX;

     step= (*pNavigatorIter)->ComputeStep( initialPosition, 
                                           initialDirection,
                                           proposedStepLength,
                                           safety ); 
     if( safety < minSafety ){ minSafety = safety; } 
     if( step < minStep )    { minStep= step; } 
     //  Later could reduce the proposed step to the latest minStep value ?
     // if( step == kInfinity ) { step = proposedStepLength; }
     fCurrentStepSize[num] = step; 
     fNewSafety[num]= safety; 
      // This is currently the safety from the last sub-step

     if( fVerboseLevel > 2 ){
       G4cout << "G4MultiNavigator::ComputeStep : Navigator [" << num << "] -- step size " << step << " safety= " << safety << G4endl;
     }
  } 

  // fWasLimitedByGeometry= false;  // <----- Could reset(?), but navigator leaves it as is
  // Whether any geometry limited the step
  // G4bool StepLimited = ( minStep <= proposedStepLength);
  // G4cout << "G4MultiNavigator::ComputeStep - StepLimited is " << StepLimited
  //  << " given minStep= " << minStep << " and proposed Step= " << proposedStepLength << G4endl;

  // Save safety value, related position
  fPreStepLocation=     initialPosition; 
  fMinSafety_PreStepPt= minSafety;

  fMinStep=  minStep; 

  G4double trueMinStep= minStep; 
  if( fMinStep == kInfinity ){
     trueMinStep = proposedStepLength;   //  Use this below for endpoint !!
  }
  fTrueMinStep = trueMinStep;

  if( fVerboseLevel > 1 ){
    G4ThreeVector endPosition;
    endPosition= initialPosition + trueMinStep * initialDirection ; 

    int oldPrec= G4cout.precision(8); 
    G4cout << "G4MultiNavigator::ComputeStep : "
           << " initialPosition = " << initialPosition 
           << " and endPosition = " << endPosition<< G4endl;
    G4cout.precision( oldPrec );
  }

  pNewSafety= minSafety; 
  // Set the EndState
  // fEndState= initialState; 
  // fEndState.SetPosition( endPosition ); 
  // fEndState.SetProperTimeOfFlight( -1.000 );   // Not defined YET
  // fEndState.SetMomentum( initialState.GetMomentum ); 

  this->WhichLimited(); 

  if( fVerboseLevel > 2 ){
    G4cout << " G4MultiNavigator::ComputeStep : exits returning " << minStep << G4endl;
  }

  return minStep;  // must return kInfinity if do not limit step
}

G4double 
G4MultiNavigator::ObtainFinalStep( G4int        navigatorId, 
                                   G4double     &pNewSafety,     // for this geom 
                                   G4double     &minStep,
                                   ELimited     &limitedStep) 
{
  G4int navigatorNo=-1; 

  if( navigatorId <= fNoActiveNavigators ){
     navigatorNo= navigatorId;
  } else { 
     G4cerr << " Navigator Id = " << navigatorId 
            << " No Active = " << fNoActiveNavigators << " . " << G4endl;
     G4Exception( "G4MultiNavigator::ObtainFinalStep : Bad Navigator Id" ); 
  }
  // if( ! ){ G4Exception( "G4MultiNavigator::ObtainFinalStep Called without call to ComputeStep"); }

  // Prepare the information to return
  pNewSafety  = fNewSafety[ navigatorNo ]; 
  limitedStep = fLimitedStep[ navigatorNo ];
  minStep= fMinStep; 

  // if( (minStep==kInfinity) || (fVerboseLevel > 1) ){ 
  if( fVerboseLevel > 1 ){ 
     G4cout << " G4MultiNavigator::ComputeStep returns " << fCurrentStepSize[ navigatorNo ]
            << " for Navigator " << navigatorNo << " Limited step = " << limitedStep 
            << " Safety(mm) = " << pNewSafety / mm << G4endl; 
  }

  return fCurrentStepSize[ navigatorNo ];
}

// ----------------------------------------------------------------------

void
G4MultiNavigator::PrepareNewTrack( const G4ThreeVector position, 
                                   const G4ThreeVector direction )
{
  if( fVerboseLevel > 1 ) 
    G4cout << " Entered G4MultiNavigator::PrepareNewTrack() " << G4endl;

  G4MultiNavigator::PrepareNavigators(); 
  //***********************************

  if( fVerboseLevel > 1 ) {
    G4cout << " Calling MultiNavigator::Locate() from G4MultiNavigator::PrepareNewTrack() " 
           << G4endl;
  }

  this->LocateGlobalPointAndSetup( position, &direction, false, false );   
  //    =========================
  // The first location for each Navigator must be non-relative
  //   or else call ResetStackAndState() for each Navigator
  // Use direction to get correct side of boundary (ignore dir= false)

  // fRelocatedPoint= false; 

  if( fVerboseLevel > 0 ) {
    G4cout << " G4MultiNavigator::PrepareNewTrack : exiting. " << G4endl;
  }

}

void
G4MultiNavigator::PrepareNavigators()
{
  // Key purposes:
  //   - Check and cache set of active navigators
  //   - Reset state for new track
  G4int num=0; 

  if( fVerboseLevel > 1 ) 
    G4cout << " G4MultiNavigator::PrepareNavigators - entered " << G4endl;
  // static G4TransportationManager* pTransportManager= 
  //       G4TransportationManager::GetTransportationManager();

  // fNavigators= true; 
  // this->MovePoint();   // Signal further that the last status is wiped

  // Message the G4NavigatorPanel / Dispatcher to find active navigators
  std::vector<G4Navigator*>::iterator pNavigatorIter; 

  fNoActiveNavigators=  pTransportManager-> GetNoActiveNavigators();
  if( fNoActiveNavigators > fMaxNav ){
    G4cerr << "Too many active Navigators (worlds). G4MultiNavigator fails." 
           << G4endl;
    G4cout << " Fatal error: Transportation Manager reports " << fNoActiveNavigators 
           << " which is more than the number allowed = "     << fMaxNav << G4endl;
    G4Exception("G4MultiNavigator::PrepareNavigators()", "TooManyNavigators",  
                FatalException,  "Too many active Navigators / worlds"); 
  }

  pNavigatorIter= pTransportManager-> GetActiveNavigatorsIterator();
  for( num=0; num< fNoActiveNavigators; ++pNavigatorIter,++num ) {
 
     // Keep information in carray ... for returning information stored
     fpNavigator[num] =  *pNavigatorIter;   
     fLimitTruth[num] = false;
     fLimitedStep[num] = kDoNot;
     fCurrentStepSize[num] = 0.0; 
     fLocatedVolume[num] = 0; 
  }
  fWasLimitedByGeometry= false; 

  // Check the world volume of the mass navigator (in case a SetWorldVolume changed it)
  G4VPhysicalVolume* massWorld =  this-> GetWorldVolume();
    // fpNavigator[0] -> GetWorldVolume();  
  if( (massWorld != fLastMassWorld) && (massWorld!=0) ) { 
     // Pass along change to Mass Navigator
     fpNavigator[0] -> SetWorldVolume( massWorld ); 
     if( fVerboseLevel > 0 ) { 
       G4cout << "G4MultiNavigator::PrepareNavigators changed world volume " 
              << " for mass geometry to " << massWorld->GetName() << G4endl; 
     }
     fLastMassWorld= massWorld;
  }else{
     if( fVerboseLevel > 2 ) { 
          G4cout << "G4MultiNavigator::PrepareNavigators retained world volume " 
                 << " Pointer= " << massWorld << G4endl;
          if( massWorld ) 
             G4cout << " Name= " << massWorld->GetName() << G4endl;
     }
  }

  if( fVerboseLevel > 2 ) {
    G4cout << " G4MultiNavigator::PrepareNavigators : exiting. " << G4endl;
  }
}


G4VPhysicalVolume* 
G4MultiNavigator::LocateGlobalPointAndSetup(const G4ThreeVector& position,
                                            const G4ThreeVector* pDirection,
                                            const G4bool pRelativeSearch,
                                            const G4bool ignoreDirection )
{
  // Locate the point in each geometry

  G4ThreeVector direction(0.0, 0.0, 0.0);
  G4bool  relative= pRelativeSearch; 
  std::vector<G4Navigator*>::iterator pNavIter= pTransportManager->GetActiveNavigatorsIterator(); 
  G4int num=0; 

  if( pDirection ) direction = *pDirection; 

#if 0
  G4ThreeVector lastEndPosition= fEndState.GetPosition(); 
  G4ThreeVector moveVec = (position - lastEndPosition );
  G4double      moveLenSq= moveVec.mag2();
  if( (!fNewTrack) && (!fRelocatedPoint) && ( moveLenSq> 0.0) ){
     ReportMove( position, lastEndPosition, "Position" ); 
     G4Exception( "G4MultiNavigator::LocateGlobalPointAndSetup", 
                  "211-LocateUnexpectedPoint",  
                  JustWarning,   
                  // FatalException,  
                  "Location is not where last ComputeStep ended."); 
  }
  fLastLocatedPosition= position; 
#endif

  if( fVerboseLevel > 2 ){
    G4cout << " G4MultiNavigator::LocateGlobalPointAndSetup : entered " << " ---------------" <<  G4endl;
    G4cout << "   Locating at position " << position << "  with direction " << direction 
           << "  relative= " << relative << " ignore direction= " << ignoreDirection<< G4endl;
    G4cout << "   Number of active navigators= " << fNoActiveNavigators << G4endl;
  }

  for ( num=0; num< fNoActiveNavigators ; ++pNavIter,++num ) {
     //  ... who limited the step ....

     // G4cout << " -- Navigator id= " << num << " NavigatorPtr " << *pNavIter << G4endl;
     // G4VPhysicalVolume* world= (*pNavIter)->GetWorldVolume(); 
     // if( world ) { G4cout << "    Navigator world= " <<  world->GetName() << G4endl; }
     // else{  G4cout << "    No world set in Navigator. " << G4endl;   }

     if( fWasLimitedByGeometry && fLimitTruth[num] ) { 
        (*pNavIter)->SetGeometricallyLimitedStep(); 
     }

     G4VPhysicalVolume *pLocated= 
     (*pNavIter)->LocateGlobalPointAndSetup( position, &direction,
     //*************************************//
                                             relative,  
                                             ignoreDirection);   
     // Set the state related to the location
     fLocatedVolume[num] = pLocated; 

     // Clear state related to the step
     fLimitedStep[num]   = kDoNot; 
     fCurrentStepSize[num] = 0.0;      
     fLimitTruth[ num ] = false;   // Always clear on locating (see Navigator)
    
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
  fWasLimitedByGeometry= false;   // Clear on locating

  if( fVerboseLevel > 2 ){
    G4cout << " G4MultiNavigator::Locate : exiting. " << G4endl << G4endl; 
  }
  // fRelocatedPoint= false;

  G4VPhysicalVolume* volMassLocated= fLocatedVolume[0]; 
  return volMassLocated;
}

void
G4MultiNavigator::LocateGlobalPointWithinVolume(const G4ThreeVector& position)
{
  // Relocate the point in each geometry
  std::vector<G4Navigator*>::iterator pNavIter= pTransportManager->GetActiveNavigatorsIterator(); 
  // const G4double cErrorTolerance=1e-12;   
  // Maximum relative error from roundoff of arithmetic 
  G4int num=0; 

  if( fVerboseLevel > 2 ){
    G4cout << G4endl; 
    G4cout << " G4MultiNavigator::ReLocate : entered " << G4endl;
    G4cout << " ----------------------   -------" <<  G4endl;
    G4cout << "  *Re*Locating at position " << position  << G4endl; 
  }

  for ( num=0; num< fNoActiveNavigators ; ++pNavIter,++num ) {
     //  ... none limited the step

     // G4VPhysicalVolume physVolume= 
     (*pNavIter)->LocateGlobalPointWithinVolume( position ); 
     //*************************************//

     // Clear state related to the step
     fLimitedStep[num]   = kDoNot; 
     fCurrentStepSize[num] = 0.0;      

     fLimitTruth[ num ] = false;   // Always clear on locating (see Navigator)
     // fLocatedVolume[num]= physVolume; 
    
     // G4cout << " ReLocated in world " << num << " at " << position << G4endl;
  }
  fWasLimitedByGeometry= false;   // Clear on locating

  fLastLocatedPosition= position; 
  // fRelocatedPoint= false;

  if( fVerboseLevel > 2 ){
    G4cout << " G4MultiNavigator::LocateGlobalPointWithinVolume : exiting " 
           << "  at position " << position << G4endl;
    G4cout << G4endl;
  }

}

// -----------------------------------------------------------------------------

G4double  G4MultiNavigator::ComputeSafety( const G4ThreeVector& position,
                                                 G4double       maxDistance)
     // Recompute safety for the relevant point
{
    G4double minSafety= DBL_MAX; 
    // G4cout << " G4MultiNavigator::ComputeSafety - called at " << position << G4endl;
  
    std::vector<G4Navigator*>::iterator pNavigatorIter;
    pNavigatorIter= pTransportManager-> GetActiveNavigatorsIterator();

    G4int num=0; 
    for( num=0; num< fNoActiveNavigators; ++pNavigatorIter,++num ) {

       G4double safety;
       safety= (*pNavigatorIter)->ComputeSafety( position, maxDistance ); 
  
       if( safety < minSafety ){ minSafety = safety; } 
       // fNewSafety[num]= safety; 
    } 

    fSafetyLocation= position;
    fMinSafety_atSafLocation = minSafety;

    if( fVerboseLevel > 1 ) { 
      G4cout << " G4MultiNavigator::ComputeSafety - returns " 
             << minSafety << " at location " << position 
             << G4endl;
    }
    return minSafety; 
}


// -----------------------------------------------------------------------------

G4TouchableHistoryHandle 
G4MultiNavigator::CreateTouchableHistoryHandle() const
{
  G4Exception( "G4MultiNavigator::CreateTouchableHistoryHandle", 
               "215-TouchableFromWrongNavigator",  
               FatalException,  
               "Getting a touchable from G4MultiNavigator is not defined."); 

  if( fVerboseLevel > 2 ){
    G4cout << "G4MultiNavigator::CreateTouchableHandle : navId = " << 0 ;
      //   << " -- " << GetNavigator(navId) << G4endl;
  }

  G4TouchableHistory* touchHist;
  touchHist= fpNavigator[0] -> CreateTouchableHistory(); 

  // G4TouchableHistory* touchHist= new G4TouchableHistory(); 

  G4VPhysicalVolume* locatedVolume= fLocatedVolume[0]; 
  if( locatedVolume == 0 )
     {
       // Workaround to ensure that the touchable is fixed !! // TODO: fix
       touchHist->UpdateYourself( locatedVolume, 
                                  touchHist->GetHistory() );
     }
    
  return G4TouchableHistoryHandle(touchHist); 
}

void
G4MultiNavigator::WhichLimited()       // Flag which processes limited the step
{
  G4int num=-1, last=-1; 
  const G4int  IdTransport= 0;  // Id of Mass Navigator !!
  G4int noLimited=0; 
  ELimited shared= kSharedOther; 

  if( fVerboseLevel > 2 )
    G4cout << " G4MultiNavigator::WhichLimited - entered " << G4endl;

  // Assume that [IdTransport] is Mass / Transport
  // G4bool transportLimited = (fCurrentStepSize[IdTransport] == fMinStep); 
  G4bool transportLimited = (fCurrentStepSize[IdTransport] == fMinStep)
                           && ( fMinStep!= kInfinity) ; 
  if( transportLimited ){ 
     shared= kSharedTransport;
  }

  for ( num= 0; num < fNoActiveNavigators; num++ ) { 
    G4bool limitedStep;

    G4double step= fCurrentStepSize[num]; 

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

#ifndef G4NO_VERBOSE
  if( fVerboseLevel > 1 ){
    this->PrintLimited();   // --> for tracing 
    G4cout << " G4MultiNavigator::WhichLimited - exiting. " << G4endl;
  }
#endif

}

void
G4MultiNavigator::PrintLimited()
{
  static G4String StrDoNot("DoNot"), StrUnique("Unique"), StrUndefined("Undefined"),
    StrSharedTransport("SharedTransport"),  StrSharedOther("SharedOther");
  // Report results -- for checking   
  G4cout << "G4MultiNavigator::PrintLimited reports: " ; 
  G4cout << "  Minimum step (true)= " << fTrueMinStep 
         << "  reported min = " << fMinStep 
         << G4endl; 
  if(  // (fCurrentStepNo <= 2) || 
      (fVerboseLevel>=2) ) {
    G4cout // << std::setw(5) << " Step#"  << " "
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
    G4cout // << std::setw(5) << fCurrentStepNo  << " " 
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

    G4Navigator *pNav= fpNavigator[ num ]; 
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
    G4cout << " G4MultiNavigator::PrintLimited - exiting. " << G4endl;
}
 

void 
G4MultiNavigator::ResetState()
{
   G4int num; 
   fWasLimitedByGeometry= false; 

   G4Exception( "G4MultiNavigator::ResetState", 
               "217-CannotImplement",  
               FatalException,  
               "Cannot call ResetState for active navigators of G4MultiNavigator.");   
   std::vector<G4Navigator*>::iterator pNavigatorIter;
   pNavigatorIter= pTransportManager-> GetActiveNavigatorsIterator();
   for( num=0; num< fNoActiveNavigators; ++pNavigatorIter,++num ) {
       //  (*pNavigatorIter)->ResetState();  // KEEP THIS comment !!!
   } 
}

void 
G4MultiNavigator::SetupHierarchy()
{
   // G4Navigator::SetupHierarchy(); 
   G4Exception( "G4MultiNavigator::SetupHierarchy", 
               "217-CannotImplement",  
               FatalException,  
               "Cannot call SetupHierarchy for active navigators of G4MultiNavigator."); 
}

void 
G4MultiNavigator::CheckMassWorld()
{
   // 
   G4VPhysicalVolume* navTrackWorld= pTransportManager->GetNavigatorForTracking()
                                                      ->GetWorldVolume();
   if( navTrackWorld != fLastMassWorld ) { 
      G4Exception( "G4MultiNavigator::CheckMassWorld", "MultiNav-220", FatalException, 
                   "Mass world pointer has been changed." ); 
   } 
}

G4VPhysicalVolume* G4MultiNavigator::ResetHierarchyAndLocate(const G4ThreeVector &point,
                                           const G4ThreeVector &direction,
                                           const G4TouchableHistory &MassHistory)
   // Reset geometry for all -- and use the touchable for the mass history
{
   G4VPhysicalVolume* massVolume=0; 
   G4int num; 
   G4Navigator* pMassNavigator= fpNavigator[0]; 

   if( pMassNavigator ){
      massVolume= pMassNavigator->ResetHierarchyAndLocate( point, direction, MassHistory); 
   }else{
      G4Exception("G4MultiNavigator::ResetHierarchyAndLocate",
                  "218-TooEarlyToReset",
                  FatalException,
                  "Cannot reset hierarchy before object is initialised with valid navigators, including a mass Navigator" );
   }

   std::vector<G4Navigator*>::iterator pNavIter= 
       pTransportManager->GetActiveNavigatorsIterator(); 

   for ( num=0; num< fNoActiveNavigators ; ++pNavIter,++num ) {
      G4bool relativeSearch, ignoreDirection; 

      (*pNavIter)-> LocateGlobalPointAndSetup( point, 
                                               &direction, 
                                               relativeSearch=false,
                                               ignoreDirection=false);
   }
   return massVolume; 
}
