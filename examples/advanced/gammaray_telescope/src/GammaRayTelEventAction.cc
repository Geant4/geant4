// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: GammaRayTelEventAction.cc,v 1.4 2000-12-06 16:53:14 flongo Exp $
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
#include "GammaRayTelPayloadHit.hh"
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
  :drawFlag("all"),trackerCollID(-1),analysisManager(aMgr)
{
}
#else
GammaRayTelEventAction::GammaRayTelEventAction()
  :drawFlag("all"), trackerCollID(-1)
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
      // This is a cycle on all the hits of this event
      for (int i=0;i<n_hit;i++)
	{
	  // Here we put the hit data in a an ASCII file for 
	  // later analysis 
	  ESil = (*CHC)[i]->GetEdepSil();
	  NStrip = (*CHC)[i]->GetNStrip();
	  NPlane = (*CHC)[i]->GetNSilPlane();
	  IsX = (*CHC)[i]->GetPlaneType();

	  outFile << G4std::setw(7) << event_id << " " << 
	    ESil/keV << " " << NStrip << 
	    " " << NPlane << " " << IsX << " " <<
	    (*CHC)[i]->GetPos().x()/mm <<" "<<
	    (*CHC)[i]->GetPos().y()/mm <<" "<<
	    (*CHC)[i]->GetPos().z()/mm <<" "<<
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













