//event action

#include "GammaRayTelEventAction.hh"
#include "GammaRayTelPayloadHit.hh"
#include "g4rw/tvordvec.h"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

// This file is a global variabe in which I store energy deposition per hit
// It is a sort of hack, I think I have to use some other way to export 
// hit information ....
extern ofstream outFile;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelEventAction::GammaRayTelEventAction()
  :drawFlag("all"), printModulo(10000),trackerCollID(-1)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelEventAction::~GammaRayTelEventAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelEventAction::BeginOfEventAction(const G4Event* evt)
{
  
 G4int evtNb = evt->GetEventID();
 if (evtNb%printModulo == 0)
   { 
    G4cout << "\n---> Event: " << evtNb << G4endl;
   }
    
 if (trackerCollID==-1)
 {
  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  trackerCollID = SDman->GetCollectionID("PayloadCollection");
 } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelEventAction::EndOfEventAction(const G4Event* evt)
{
  G4int event_id = evt->GetEventID();

  G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();


  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  GammaRayTelPayloadHitsCollection* CHC = NULL;
  if (HCE)
      CHC = (GammaRayTelPayloadHitsCollection*)(HCE->GetHC(trackerCollID));
       
    if (CHC)
      {
	int n_hit = CHC->entries();
	G4cout << "Number of hits in this event =  " << n_hit << G4endl;
	G4double ESil=0;
	G4int NStrip, NPlane, IsX;
	for (int i=0;i<n_hit;i++)
	  {
	    // Here I put the energy deposition per hit in the file for 
	    // later analysis 
	    ESil = (*CHC)[i]->GetEdepSil();
	    NStrip = (*CHC)[i]->GetNStrip();
	    NPlane = (*CHC)[i]->GetNSilPlane();
	    IsX = (*CHC)[i]->GetPlaneType();
	    outFile << G4std::setw(7) << ESil/keV << " " << NStrip << 
		    " " << NPlane << " " << IsX << G4endl;
	  }
      }
    
    
    
    if(G4VVisManager::GetConcreteInstance())
  {
    for(G4int i=0; i<n_trajectories; i++) 
         { G4Trajectory* trj = (G4Trajectory *)((*(evt->GetTrajectoryContainer()))[i]);
           if (drawFlag == "all") trj->DrawTrajectory(50);
           else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.))
                                  trj->DrawTrajectory(50); 
         }
  }
}










