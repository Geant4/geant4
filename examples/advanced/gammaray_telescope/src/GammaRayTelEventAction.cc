// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: GammaRayTelEventAction.cc,v 1.5 2001-03-05 13:58:23 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ GammaRayTelEventAction  ------
//           by  R.Giannitrapani, F.Longo & G.Santin (13 nov 2000)
//
// ************************************************************

#include "GammaRayTelEventAction.hh"
#include "GammaRayTelTrackerHit.hh"
#include "GammaRayTelAnticoincidenceHit.hh"
#include "GammaRayTelCalorimeterHit.hh"

#include "g4rw/tvordvec.h"

#ifdef G4ANALYSIS_USE
#include "GammaRayTelAnalysisManager.hh"
#endif

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

// This file is a global variable in which we store energy deposition per hit
// and other relevant information
extern G4std::ofstream outFile;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifdef G4ANALYSIS_USE
GammaRayTelEventAction::GammaRayTelEventAction(GammaRayTelAnalysisManager* aMgr)
  :drawFlag("all"),trackerCollID(-1),  calorimeterCollID(-1),                
  anticoincidenceCollID(-1), analysisManager(aMgr)
{
}
#else
GammaRayTelEventAction::GammaRayTelEventAction()
  :drawFlag("all"), trackerCollID(-1),calorimeterCollID(-1),                
  anticoincidenceCollID(-1)
{
}
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelEventAction::~GammaRayTelEventAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelEventAction::BeginOfEventAction(const G4Event* evt)
{

  G4int evtNb = evt->GetEventID();
  G4cout << "Event: " << evtNb << G4endl;

  G4SDManager * SDman = G4SDManager::GetSDMpointer();  
  if (trackerCollID==-1)
    {
      trackerCollID = SDman->GetCollectionID("TrackerCollection");
    }
  if(anticoincidenceCollID==-1)
    {
      anticoincidenceCollID =
	SDman->GetCollectionID("AnticoincidenceCollection");
    }
  if(calorimeterCollID==-1)
    {
      calorimeterCollID =
	SDman->GetCollectionID("CalorimeterCollection");
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
  GammaRayTelTrackerHitsCollection* THC = 0;
  GammaRayTelCalorimeterHitsCollection* CHC = 0;
  GammaRayTelAnticoincidenceHitsCollection* AHC = 0;


  if (HCE)
    {
      THC = (GammaRayTelTrackerHitsCollection*)(HCE->GetHC(trackerCollID));
      CHC = (GammaRayTelCalorimeterHitsCollection*)
	(HCE->GetHC(calorimeterCollID));
      AHC = (GammaRayTelAnticoincidenceHitsCollection*)
	(HCE->GetHC(anticoincidenceCollID));
    }
  
  if (THC)
    {
      int n_hit = THC->entries();
      G4cout << "Number of tracker hits in this event =  " << n_hit << G4endl;
      G4double ESil=0;
      G4int NStrip, NPlane, IsX;
      
      // This is a cycle on all the tracker hits of this event

      for (int i=0;i<n_hit;i++)
	{
	  // Here we put the hit data in a an ASCII file for 
	  // later analysis 
	  ESil = (*THC)[i]->GetEdepSil();
	  NStrip = (*THC)[i]->GetNStrip();
	  NPlane = (*THC)[i]->GetNSilPlane();
	  IsX = (*THC)[i]->GetPlaneType();
	  
	  outFile << G4std::setw(7) << event_id << " " << 
	    ESil/keV << " " << NStrip << 
	    " " << NPlane << " " << IsX << " " <<
	    (*THC)[i]->GetPos().x()/mm <<" "<<
	    (*THC)[i]->GetPos().y()/mm <<" "<<
	    (*THC)[i]->GetPos().z()/mm <<" "<<
	    G4endl;
	  
#ifdef G4ANALYSIS_USE
	  // Here we fill the histograms of the Analysis manager
	  if(IsX)
	    {
	      if (analysisManager->GetHisto2DMode()=="position")
		analysisManager->InsertPositionXZ((*CHC)[i]->GetPos().x()/mm,(*CHC)[i]->GetPos().z()/mm);
	      else
		analysisManager->InsertPositionXZ(NStrip, NPlane);  	      
	      if (NPlane == 0) analysisManager->InsertEnergy(ESil/keV);
	      analysisManager->InsertHits(NPlane);
	    }
	  else
	    if (analysisManager->GetHisto2DMode()=="position")
	      analysisManager->InsertPositionYZ((*CHC)[i]->GetPos().y()/mm,(*CHC)[i]->GetPos().z()/mm);  
	    else 
	      analysisManager->InsertPositionYZ(NStrip, NPlane);  	      
#endif
	}

      // Here we call the analysis manager function for visualization
#ifdef G4ANALYSIS_USE
      analysisManager->EndOfEvent(n_hit);
#endif
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













