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
// $Id: GammaRayTelEventAction.cc,v 1.15 2002-11-14 10:55:28 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTelEventAction  ------
//           by  R.Giannitrapani, F.Longo & G.Santin (13 nov 2000)
//
// - inclusion of Digits by F.Longo & R.Giannitrapani (24 oct 2001)
// 
// - Modification of analysis management by G.Santin (18 Nov 2001)
// 
// ************************************************************

#include "GammaRayTelEventAction.hh"
#include "GammaRayTelTrackerHit.hh"
#include "GammaRayTelAnticoincidenceHit.hh"
#include "GammaRayTelCalorimeterHit.hh"

#ifdef G4ANALYSIS_USE
#include "GammaRayTelAnalysis.hh"
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

#include "GammaRayTelDigi.hh"
#include "GammaRayTelDigitizer.hh"
#include "G4DigiManager.hh"

// This file is a global variable in which we store energy deposition per hit
// and other relevant information

#ifdef G4STORE_DATA
extern G4std::ofstream outFile;
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelEventAction::GammaRayTelEventAction()
  :trackerCollID(-1),calorimeterCollID(-1),                
  anticoincidenceCollID(-1), drawFlag("all")
{ 
  G4DigiManager * fDM = G4DigiManager::GetDMpointer();
  GammaRayTelDigitizer * myDM = new GammaRayTelDigitizer( "TrackerDigitizer" );
  fDM->AddNewModule(myDM);
}

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

  if (trackerCollID==-1) {
    trackerCollID = SDman->GetCollectionID("TrackerCollection");
  }
  if(anticoincidenceCollID==-1) {
    anticoincidenceCollID =
      SDman->GetCollectionID("AnticoincidenceCollection");
  }
  if(calorimeterCollID==-1) {
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


  G4DigiManager * fDM = G4DigiManager::GetDMpointer();

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
	  
#ifdef G4STORE_DATA
	  outFile << G4std::setw(7) << event_id << " " << 
	    ESil/keV << " " << NStrip << 
	    " " << NPlane << " " << IsX << " " <<
	    (*THC)[i]->GetPos().x()/mm <<" "<<
	    (*THC)[i]->GetPos().y()/mm <<" "<<
	    (*THC)[i]->GetPos().z()/mm <<" "<<
	    G4endl;
#else 	  
	  G4cout << G4std::setw(7) << event_id << " " << 
	    ESil/keV << " " << NStrip << 
	    " " << NPlane << " " << IsX << " " <<
	    (*THC)[i]->GetPos().x()/mm <<" "<<
	    (*THC)[i]->GetPos().y()/mm <<" "<<
	    (*THC)[i]->GetPos().z()/mm <<" "<<
	    G4endl;
#endif

#ifdef G4ANALYSIS_USE

	  // Here we fill the histograms of the Analysis manager
	  GammaRayTelAnalysis* analysis = GammaRayTelAnalysis::getInstance();

	  if(IsX) 
	    {
	      if (analysis->GetHisto2DMode()=="position")
		analysis->InsertPositionXZ((*THC)[i]->GetPos().x()/mm,(*THC)[i]->GetPos().z()/mm);
	      else
		analysis->InsertPositionXZ(NStrip, NPlane);  	      
	      if (NPlane == 0) analysis->InsertEnergy(ESil/keV);
	      analysis->InsertHits(NPlane);
	    } 
	  else 
	    {
	      if (analysis->GetHisto2DMode()=="position")
		analysis->InsertPositionYZ((*THC)[i]->GetPos().y()/mm,(*THC)[i]->GetPos().z()/mm);  
	      else 
		analysis->InsertPositionYZ(NStrip, NPlane);  	      
	      if (NPlane == 0) analysis->InsertEnergy(ESil/keV);
	      analysis->InsertHits(NPlane);
	    }
	  
#ifdef G4ANALYSIS_USE_NTUPLE
	  analysis->setNtuple( ESil/keV, NPlane, (*THC)[i]->GetPos().x()/mm,
			       (*THC)[i]->GetPos().y()/mm,
			       (*THC)[i]->GetPos().z()/mm);
#endif
	  
#endif
	  
	}
      // Here we call the analysis manager function for visualization
#ifdef G4ANALYSIS_USE
      GammaRayTelAnalysis* analysis = GammaRayTelAnalysis::getInstance();
      analysis->EndOfEvent(n_hit);
#endif
    }
  
  GammaRayTelDigitizer * myDM = 
    (GammaRayTelDigitizer*)fDM->FindDigitizerModule( "TrackerDigitizer" );
  myDM->Digitize();
  
  G4int myDigiCollID = fDM->GetDigiCollectionID("DigitsCollection");

  // G4cout << "digi collecion" << myDigiCollID << G4endl;
  
  GammaRayTelDigitsCollection * DC = (GammaRayTelDigitsCollection*)fDM->GetDigiCollection( myDigiCollID );
  
  if(DC) {
    //    G4cout << "Total Digits " << DC->entries() << G4endl;
    G4int n_digi =  DC->entries();
    G4int NStrip, NPlane, IsX;
    for (G4int i=0;i<n_digi;i++) {
      // Here we put the digi data in a an ASCII file for 
      // later analysis
      NStrip = (*DC)[i]->GetStripNumber();
      NPlane = (*DC)[i]->GetPlaneNumber();
      IsX = (*DC)[i]->GetPlaneType();
      
      //      outFile << G4std::setw(7) << event_id << " " << NStrip << 
      //	" " << NPlane << " " << IsX << " " << G4endl;	
    }
  }
  
  if(G4VVisManager::GetConcreteInstance()) {
    for(G4int i=0; i<n_trajectories; i++) { 
      G4Trajectory* trj = (G4Trajectory *)((*(evt->GetTrajectoryContainer()))[i]);
      if (drawFlag == "all") trj->DrawTrajectory(50);
      else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.))
	trj->DrawTrajectory(50); 
    }
  }
}













