// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: NTSTEventAction.cc,v 1.2 2003-11-07 22:08:59 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4Timer.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4VTrajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "NTSTEventAction.hh"
#include "NTSTEventActionMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

NTSTEventAction::NTSTEventAction()
  : EventTime(new G4Timer()),
   MeanUserEventTime(0),  RmsUserEventTime(0),
   MeanRealEventTime(0),  RmsRealEventTime(0), 
   NumberOfEvents(0), 
   MeanVertices(0), RmsVertices(0),
   MeanTracks(0), RmsTracks(0), 
   drawFlag("all"), eventMessenger(NULL)
{
  eventMessenger = new NTSTEventActionMessenger(this);
}

#include <iomanip.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

NTSTEventAction::~NTSTEventAction()
{
  if (NumberOfEvents>0) {
    G4cout << "### Processed number of events for all runs: " 
	   << NumberOfEvents << G4endl;
    MeanUserEventTime = MeanUserEventTime / NumberOfEvents;
    MeanRealEventTime = MeanRealEventTime / NumberOfEvents;
    RmsUserEventTime = RmsUserEventTime / NumberOfEvents;
    RmsRealEventTime = RmsRealEventTime / NumberOfEvents;
    G4double ErrUserEventTime = 
      sqrt((RmsUserEventTime - MeanUserEventTime*MeanUserEventTime)
	   /NumberOfEvents);
    G4double ErrRealEventTime = 
      sqrt((RmsRealEventTime - MeanRealEventTime*MeanRealEventTime)
	   /NumberOfEvents);
    G4cout << setprecision(3) 
	   << "### Event user time = " << setw(6) << MeanUserEventTime
	   << setprecision(3) 
	   << " +- " << setw(6) << ErrUserEventTime << " (sec) " << G4endl;   
    G4cout << setprecision(3)
	   << "### Event real time = " << setw(6) << MeanRealEventTime
	   << setprecision(3) 
	   << " +- " << setw(6) << ErrRealEventTime << " (sec) " << G4endl;   
  
    MeanVertices = MeanVertices / NumberOfEvents;
    RmsVertices = RmsVertices / NumberOfEvents;
    MeanTracks = MeanTracks / NumberOfEvents;
    RmsTracks = RmsTracks / NumberOfEvents;
    G4double ErrVertices = 
      sqrt((RmsVertices - MeanVertices*MeanVertices)/NumberOfEvents);
    G4double ErrTracks =
      sqrt((RmsTracks - MeanTracks*MeanTracks)/NumberOfEvents);

    G4cout << setprecision(3)
	   << "### Number of Vertices = " << setw(6) << MeanVertices << " +- "
	   << setw(6) << ErrVertices
	   << " Number of Tracks = " << setw(6) << MeanTracks << " +- "
	   << setw(6) << ErrTracks << G4endl;

				
  }
  delete eventMessenger;
  delete EventTime;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void NTSTEventAction::BeginOfEventAction(const G4Event* ) // evt)
{  
  EventTime->Start();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void NTSTEventAction::EndOfEventAction(const G4Event* evt)
{
  EventTime->Stop();
  G4cout << "### Event " << setprecision(3) 
	 << evt->GetEventID()+1 << " " << *EventTime << G4endl;

  // event statistics
  G4double ElapsedUserTime = EventTime->GetUserElapsed();
  G4double ElapsedRealTime = EventTime->GetRealElapsed();
  MeanUserEventTime += ElapsedUserTime;
  MeanRealEventTime += ElapsedRealTime;
  RmsUserEventTime += ElapsedUserTime*ElapsedUserTime;
  RmsRealEventTime += ElapsedRealTime*ElapsedRealTime;
  NumberOfEvents++;

  // vertex, track statistics
  G4int Vertices = evt->GetNumberOfPrimaryVertex();
  MeanVertices+=Vertices;
  RmsVertices+=Vertices*Vertices;
  G4int Tracks=0;
  for (G4int iv=0; iv<Vertices; iv++){
    Tracks+=evt->GetPrimaryVertex(iv)->GetNumberOfParticle();
  }
  MeanTracks+=Tracks;
  RmsTracks+=Tracks*Tracks;

  G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if(trajectoryContainer)
    { n_trajectories = trajectoryContainer->entries(); }
  G4cout << "    " << n_trajectories 
	 << " trajectories stored in this event." << G4endl;


  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if(pVVisManager) {
    for(G4int i=0; i<n_trajectories; i++) {
      G4VTrajectory* trj = (*(evt->GetTrajectoryContainer()))[i];
      //           if (drawFlag == "all") trj->DrawTrajectory(50);
      //           else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.))
      trj->DrawTrajectory(50); 
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
