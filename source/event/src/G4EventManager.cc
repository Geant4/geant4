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
// $Id: G4EventManager.cc 104522 2017-06-02 07:19:30Z gcosmo $
//
//
//

#include "G4EventManager.hh"
#include "G4ios.hh"
#include "G4EvManMessenger.hh"
#include "G4Event.hh"
#include "G4UserEventAction.hh"
#include "G4UserStackingAction.hh"
#include "G4SDManager.hh"
#include "G4StateManager.hh"
#include "G4ApplicationState.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "Randomize.hh"

G4ThreadLocal G4EventManager* G4EventManager::fpEventManager = nullptr;
G4EventManager* G4EventManager::GetEventManager()
{ return fpEventManager; }

G4EventManager::G4EventManager()
:currentEvent(nullptr),trajectoryContainer(nullptr),
 trackIDCounter(0),
 verboseLevel(0),tracking(false),abortRequested(false),
 storetRandomNumberStatusToG4Event(false)
{
 if(fpEventManager)
 {
  G4Exception("G4EventManager::G4EventManager","Event0001",FatalException,
  "G4EventManager::G4EventManager() has already been made.");
 }
 else
 {
  trackManager = new G4TrackingManager;
  transformer = new G4PrimaryTransformer;
  trackContainer = new G4StackManager;
  theMessenger = new G4EvManMessenger(this);
  sdManager = G4SDManager::GetSDMpointerIfExist();
  stateManager = G4StateManager::GetStateManager();
  fpEventManager = this;
  userEventAction = nullptr;
  userStackingAction = nullptr;
  userTrackingAction = nullptr;
  userSteppingAction = nullptr;
 }
}

/* private -> never called
G4EventManager::G4EventManager(const G4EventManager&) {;}
G4EventManager& G4EventManager::operator=(const G4EventManager&)
{ return *this; }
*/

G4EventManager::~G4EventManager()
{
   delete trackContainer;
   delete transformer;
   delete trackManager;
   delete theMessenger;
   delete userEventAction;
   fpEventManager = 0;
}

/*
const G4EventManager & G4EventManager::operator=(const G4EventManager &right)
{ }
G4int G4EventManager::operator==(const G4EventManager &right) const { }
G4int G4EventManager::operator!=(const G4EventManager &right) const { }
*/

void G4EventManager::DoProcessing(G4Event* anEvent)
{
  abortRequested = false;
  G4ApplicationState currentState = stateManager->GetCurrentState();
  if(currentState!=G4State_GeomClosed)
  {
    G4Exception("G4EventManager::ProcessOneEvent",
                "Event0002", JustWarning,
                "IllegalApplicationState -- Geometry is not closed : cannot process an event.");
    return;
  }
  currentEvent = anEvent;
  stateManager->SetNewState(G4State_EventProc);
  if(storetRandomNumberStatusToG4Event>1)
  {
    std::ostringstream oss;
    CLHEP::HepRandom::saveFullState(oss);
    randomNumberStatusToG4Event = oss.str();
    currentEvent->SetRandomNumberStatusForProcessing(randomNumberStatusToG4Event); 
  }

  // Reseting Navigator has been moved to G4EventManager, so that resetting
  // is now done for every event.
  G4ThreeVector center(0,0,0);
  G4Navigator* navigator =
      G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  navigator->LocateGlobalPointAndSetup(center,0,false);
                                                                                      
  G4Track * track = nullptr;
  G4TrackStatus istop = fAlive;

#ifdef G4VERBOSE
  if ( verboseLevel > 0 )
  {
    G4cout << "=====================================" << G4endl;
    G4cout << "  G4EventManager::ProcessOneEvent()  " << G4endl;
    G4cout << "=====================================" << G4endl;
  }
#endif

  trackContainer->PrepareNewEvent();

#ifdef G4_STORE_TRAJECTORY
  trajectoryContainer = nullptr;
#endif

  sdManager = G4SDManager::GetSDMpointerIfExist();
  if(sdManager)
  { currentEvent->SetHCofThisEvent(sdManager->PrepareNewEvent()); }

  if(userEventAction) userEventAction->BeginOfEventAction(currentEvent);

#ifdef G4VERBOSE
  if ( verboseLevel > 1 )
  {
    G4cout << currentEvent->GetNumberOfPrimaryVertex()
         << " vertices passed from G4Event." << G4endl;
  }
#endif

  if(!abortRequested)
  { StackTracks( transformer->GimmePrimaries( currentEvent, trackIDCounter ),true ); }

#ifdef G4VERBOSE
  if ( verboseLevel > 0 )
  {
    G4cout << trackContainer->GetNTotalTrack() << " primaries "
         << "are passed from G4EventTransformer." << G4endl;
    G4cout << "!!!!!!! Now start processing an event !!!!!!!" << G4endl;
  }
#endif
  
  G4VTrajectory* previousTrajectory;
  while( ( track = trackContainer->PopNextTrack(&previousTrajectory) ) != 0 ) // Loop checking 12.28.2015 M.Asai
  {

#ifdef G4VERBOSE
    if ( verboseLevel > 1 )
    {
      G4cout << "Track " << track << " (trackID " << track->GetTrackID()
  	 << ", parentID " << track->GetParentID() 
  	 << ") is passed to G4TrackingManager." << G4endl;
    }
#endif

    tracking = true;
    trackManager->ProcessOneTrack( track );
    istop = track->GetTrackStatus();
    tracking = false;

#ifdef G4VERBOSE
    if ( verboseLevel > 0 )
    {
      G4cout << "Track (trackID " << track->GetTrackID()
	 << ", parentID " << track->GetParentID()
         << ") is processed with stopping code " << istop << G4endl;
    }
#endif

    G4VTrajectory * aTrajectory = nullptr;
#ifdef G4_STORE_TRAJECTORY
    aTrajectory = trackManager->GimmeTrajectory();

    if(previousTrajectory)
    {
      previousTrajectory->MergeTrajectory(aTrajectory);
      delete aTrajectory;
      aTrajectory = previousTrajectory;
    }
    if(aTrajectory&&(istop!=fStopButAlive)&&(istop!=fSuspend))
    {
      if(!trajectoryContainer)
      { trajectoryContainer = new G4TrajectoryContainer; 
        currentEvent->SetTrajectoryContainer(trajectoryContainer); }
      trajectoryContainer->insert(aTrajectory);
    }
#endif

    G4TrackVector * secondaries = trackManager->GimmeSecondaries();
    switch (istop)
    {
      case fStopButAlive:
      case fSuspend:
        trackContainer->PushOneTrack( track, aTrajectory );
        StackTracks( secondaries );
        break;

      case fPostponeToNextEvent:
        trackContainer->PushOneTrack( track );
        StackTracks( secondaries );
        break;

      case fStopAndKill:
        StackTracks( secondaries );
        delete track;
        break;

      case fAlive:
        G4Exception("G4EventManager::DoProcessing","Event004",JustWarning,
            "Illegal trackstatus returned from G4TrackingManager. Continue with"\
            "simulation.");
        break;
      case fKillTrackAndSecondaries:
        //if( secondaries ) secondaries->clearAndDestroy();
        if( secondaries )
        {
          for(size_t i=0;i<secondaries->size();i++)
          { delete (*secondaries)[i]; }
          secondaries->clear();
        }
        delete track;
        break;
    }
  }

#ifdef G4VERBOSE
  if ( verboseLevel > 0 )
  {
    G4cout << "NULL returned from G4StackManager." << G4endl;
    G4cout << "Terminate current event processing." << G4endl;
  }
#endif

  if(sdManager)
  { sdManager->TerminateCurrentEvent(currentEvent->GetHCofThisEvent()); }

  if(userEventAction) userEventAction->EndOfEventAction(currentEvent);

  stateManager->SetNewState(G4State_GeomClosed);
  currentEvent = nullptr;
  abortRequested = false;
}

void G4EventManager::StackTracks(G4TrackVector *trackVector,G4bool IDhasAlreadySet)
{
  if( trackVector )
  {
    //size_t n_passedTrack = trackVector->size();
    //if( n_passedTrack == 0 ) return;
    //for( size_t i = 0; i < n_passedTrack; i++ )
    //{
    //  newTrack = (*trackVector)[ i ];
    if( trackVector->size() == 0 ) return;
    for( auto newTrack : *trackVector )
    {
      trackIDCounter++;
      if(!IDhasAlreadySet)
      {
        newTrack->SetTrackID( trackIDCounter );
        if(newTrack->GetDynamicParticle()->GetPrimaryParticle())
        {
          G4PrimaryParticle* pp
            = (G4PrimaryParticle*)(newTrack->GetDynamicParticle()->GetPrimaryParticle());
          pp->SetTrackID(trackIDCounter);
        }
      }
      newTrack->SetOriginTouchableHandle(newTrack->GetTouchableHandle());
      trackContainer->PushOneTrack( newTrack );
#ifdef G4VERBOSE
      if ( verboseLevel > 1 )
      {
        G4cout << "A new track " << newTrack 
             << " (trackID " << newTrack->GetTrackID()
	     << ", parentID " << newTrack->GetParentID() 
	     << ") is passed to G4StackManager." << G4endl;
      }
#endif
    }
    trackVector->clear();
  }
}

void G4EventManager::SetUserAction(G4UserEventAction* userAction)
{
  userEventAction = userAction;
  if(userEventAction) userEventAction->SetEventManager(this);
}

void G4EventManager::SetUserAction(G4UserStackingAction* userAction)
{
  userStackingAction = userAction;
  trackContainer->SetUserStackingAction(userAction);
}

void G4EventManager::SetUserAction(G4UserTrackingAction* userAction)
{
  userTrackingAction = userAction;
  trackManager->SetUserAction(userAction);
}

void G4EventManager::SetUserAction(G4UserSteppingAction* userAction)
{
  userSteppingAction = userAction;
  trackManager->SetUserAction(userAction);
}

void G4EventManager::ProcessOneEvent(G4Event* anEvent)
{
  trackIDCounter = 0;
  DoProcessing(anEvent);
}

void G4EventManager::ProcessOneEvent(G4TrackVector* trackVector,G4Event* anEvent)
{
  static G4ThreadLocal G4String *randStat = 0;
  if (!randStat) randStat = new G4String;
  trackIDCounter = 0;
  G4bool tempEvent = false;
  if(!anEvent)
  {
    anEvent = new G4Event();
    tempEvent = true;
  }
  if(storetRandomNumberStatusToG4Event==1 || storetRandomNumberStatusToG4Event==3)
  {
    std::ostringstream oss;
    CLHEP::HepRandom::saveFullState(oss);
    anEvent->SetRandomNumberStatus(*randStat=oss.str());
  }
  StackTracks(trackVector,false);
  DoProcessing(anEvent);
  if(tempEvent)
  { delete anEvent; }
}

void G4EventManager::SetUserInformation(G4VUserEventInformation* anInfo)
{ 
  G4ApplicationState currentState = stateManager->GetCurrentState();
  if(currentState!=G4State_EventProc || currentEvent==0)
  {
    G4Exception("G4EventManager::SetUserInformation",
                "Event0003", JustWarning,
                "G4VUserEventInformation cannot be set because of ansense of G4Event.");
    return;
  }
  
  currentEvent->SetUserInformation(anInfo);
}

G4VUserEventInformation* G4EventManager::GetUserInformation()
{ 
  G4ApplicationState currentState = stateManager->GetCurrentState();
  if(currentState!=G4State_EventProc || currentEvent==0)
  { return 0; }
  
  return currentEvent->GetUserInformation();
}

void G4EventManager::KeepTheCurrentEvent()
{ if(currentEvent) currentEvent->KeepTheEvent(); }

void G4EventManager::AbortCurrentEvent()
{
  abortRequested = true;
  trackContainer->clear();
  if(tracking) trackManager->EventAborted();
}

