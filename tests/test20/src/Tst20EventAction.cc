// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst20EventAction.cc,v 1.1 2001-05-24 19:49:30 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ Tst20EventAction  ------
//           by  R.Giannitrapani, F.Longo & G.Santin (13 nov 2000)
//
// ************************************************************

#include "Tst20EventAction.hh"
#include "Tst20TrackerHit.hh"
#include "Tst20AnticoincidenceHit.hh"
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst20EventAction::Tst20EventAction()
  :drawFlag("all"), trackerCollID(-1), anticoincidenceCollID(-1)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst20EventAction::~Tst20EventAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst20EventAction::BeginOfEventAction(const G4Event* evt)
{
  G4int evtNb = evt->GetEventID();

  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  
  if(trackerCollID<0||anticoincidenceCollID<0)
    {
      trackerCollID = SDman->GetCollectionID("TrackerCollection");
      anticoincidenceCollID =
	SDman->GetCollectionID("AnticoincidenceCollection");
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst20EventAction::EndOfEventAction(const G4Event* evt)
{
  G4int event_id = evt->GetEventID();
  
  if(trackerCollID<0||anticoincidenceCollID<0) return;
  
  G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  Tst20TrackerHitsCollection* THC = NULL;
  Tst20AnticoincidenceHitsCollection* AHC = NULL;
  
  
  if (HCE)
    {
      THC = (Tst20TrackerHitsCollection*)(HCE->GetHC(trackerCollID));
      AHC = (Tst20AnticoincidenceHitsCollection*)
	(HCE->GetHC(anticoincidenceCollID));
    }
  
  if (THC)
    {
      int n_hit = THC->entries();

      G4double ESil=0;
      G4int NPixel, NPlane;
      // This is a cycle on all the hits of this event

      for (int i=0;i<n_hit;i++)
	{

	  ESil = (*THC)[i]->GetEdepSil();
	  NPixel = (*THC)[i]->GetNPixel();
	  NPlane = (*THC)[i]->GetNSilPlane();

	  G4cout << G4std::setw(7) << event_id << " " << 
	    ESil/keV << " " << NPixel << 
	    " " << NPlane << " " <<
	    (*THC)[i]->GetPos().x()/mm <<" "<<
	    (*THC)[i]->GetPos().y()/mm <<" "<<
	    (*THC)[i]->GetPos().z()/mm <<" "<<
	    G4endl;
	  
	}
      
    }

  if(AHC)
    {
      int n_hit = AHC->entries();
      G4double totE = 0;
      for(int i=0;i<n_hit;i++)
	{ 
	  totE += (*AHC)[i]->GetEdepACD(); 
	}
      if(totE>0.)
	{
	  G4cout << G4std::setw(7) << event_id << 
	    "  Total energy deposition in anticoincidence : " << 
	    totE / keV << " (keV)" << G4endl;
	}      
    }
  
}



































