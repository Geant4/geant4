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
// $Id: G4EventManager.cc,v 1.11 2001-08-28 05:35:00 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

G4EventManager* G4EventManager::fpEventManager = 0;
G4EventManager* G4EventManager::GetEventManager()
{ return fpEventManager; }

G4EventManager::G4EventManager()
:currentEvent(0),trajectoryContainer(0),
 verboseLevel(0),tracking(false)
{
 if(fpEventManager)
 {
  G4Exception("G4EventManager::G4EventManager() has already been made.");
 }
 else
 {
  trackManager = new G4TrackingManager;
  transformer = new G4PrimaryTransformer;
  trackContainer = new G4StackManager;
  theMessenger = new G4EvManMessenger(this);
  sdManager = G4SDManager::GetSDMpointerIfExist();
  fpEventManager = this;
  userEventAction = 0;
  userStackingAction = 0;
  userTrackingAction = 0;
  userSteppingAction = 0;
 }
}

// private -> never called
G4EventManager::G4EventManager(const G4EventManager &right) { }
G4EventManager& G4EventManager::operator=(const G4EventManager& right)
{ 
  return *this;
}

G4EventManager::~G4EventManager()
{
   delete trackContainer;
   delete transformer;
   delete trackManager;
   delete theMessenger;
   if(userEventAction) delete userEventAction;
   fpEventManager = 0;
}

/*
const G4EventManager & G4EventManager::operator=(const G4EventManager &right)
{ }
G4int G4EventManager::operator==(const G4EventManager &right) const { }
G4int G4EventManager::operator!=(const G4EventManager &right) const { }
*/


void G4EventManager::ProcessOneEvent(G4Event* anEvent)
{
  currentEvent = anEvent;
  G4Track * track;
  G4TrackStatus istop;

#ifdef G4VERBOSE
  if ( verboseLevel > 0 )
  {
    G4cout << "=====================================" << G4endl;
    G4cout << "  G4EventManager::ProcessOneEvent()  " << G4endl;
    G4cout << "=====================================" << G4endl;
  }
#endif

  trackIDCounter = trackContainer->PrepareNewEvent();

#ifdef G4_STORE_TRAJECTORY
  trajectoryContainer = 0;
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

  StackTracks( transformer->GimmePrimaries( currentEvent, trackIDCounter ),true );

#ifdef G4VERBOSE
  if ( verboseLevel > 0 )
  {
    G4cout << trackContainer->GetNTotalTrack() << " primaries "
         << "are passed from G4EventTransformer." << G4endl;
    G4cout << "!!!!!!! Now start processing an event !!!!!!!" << G4endl;
  }
#endif
  
  G4VTrajectory* previousTrajectory;
  while( ( track = trackContainer->PopNextTrack(&previousTrajectory) ) != 0 )
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

    G4VTrajectory * aTrajectory = 0;
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
      { trajectoryContainer = new G4TrajectoryContainer; }
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
        G4cout << "Illeagal TrackStatus returned from G4TrackingManager!"
             << G4endl;
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

#ifdef G4_STORE_TRAJECTORY
  currentEvent->SetTrajectoryContainer(trajectoryContainer);
#endif

  if(userEventAction) userEventAction->EndOfEventAction(currentEvent);
  currentEvent = 0;

}

void G4EventManager::StackTracks(G4TrackVector *trackVector,G4bool IDhasAlreadySet)
{
  G4Track * newTrack;

  if( trackVector )
  {

    size_t n_passedTrack = trackVector->size();
    if( n_passedTrack == 0 ) return;
    for( size_t i = 0; i < n_passedTrack; i++ )
    {
      newTrack = (*trackVector)[ i ];
      trackIDCounter++;
      if(!IDhasAlreadySet) newTrack->SetTrackID( trackIDCounter );
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


