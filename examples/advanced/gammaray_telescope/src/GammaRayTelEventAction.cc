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

#include "GammaRayTelAnalysis.hh"

#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelEventAction::GammaRayTelEventAction(GammaRayTelRunAction* runAction)
  :trackerCollID(-1),calorimeterCollID(-1),                
   anticoincidenceCollID(-1), drawFlag("all"),
   theRunAction(runAction)
{ 
  GammaRayTelDigitizer * myDM = new GammaRayTelDigitizer( "GammaRayTelDigitizer" );
  G4DigiManager::GetDMpointer()->AddNewModule(myDM);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelEventAction::~GammaRayTelEventAction()
{;}

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


  if (!theRunAction)
    {
      G4Exception("GammaRayTelEventAction::BeginOfEventAction()",
		  "GTR0001",FatalException,"Null pointer to Run Action: this should not be");
    }
#ifdef G4STORE_DATA
  std::ofstream* outFile = theRunAction->GetOutputFile();
#endif

  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  GammaRayTelTrackerHitsCollection* THC = 0;

  //  GammaRayTelCalorimeterHitsCollection* CHC = 0;
  //  GammaRayTelAnticoincidenceHitsCollection* AHC = 0;


  G4DigiManager * fDM = G4DigiManager::GetDMpointer();
  GammaRayTelAnalysis* analysis = GammaRayTelAnalysis::getInstance();

  if (HCE)
    {
      THC = (GammaRayTelTrackerHitsCollection*)(HCE->GetHC(trackerCollID));
      //      CHC = (GammaRayTelCalorimeterHitsCollection*)
      //	(HCE->GetHC(calorimeterCollID));
      //      AHC = (GammaRayTelAnticoincidenceHitsCollection*)
      //	(HCE->GetHC(anticoincidenceCollID));
      
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
	     

	      (*outFile) << std::setw(7) << event_id << " " << 
		ESil/keV << " " << NStrip << 
		" " << NPlane << " " << IsX << " " <<
		(*THC)[i]->GetPos().x()/mm <<" "<<
		(*THC)[i]->GetPos().y()/mm <<" "<<
		(*THC)[i]->GetPos().z()/mm <<" "<<
		G4endl;
#else 	  
	      G4cout << std::setw(7) << event_id << " " << 
		ESil/keV << " " << NStrip << 
		" " << NPlane << " " << IsX << " " <<
		(*THC)[i]->GetPos().x()/mm <<" "<<
		(*THC)[i]->GetPos().y()/mm <<" "<<
		(*THC)[i]->GetPos().z()/mm <<" "<<
		G4endl;
#endif
	      
	      
	      // Here we fill the histograms of the Analysis manager
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
	      
	      analysis->setNtuple( ESil/keV, NPlane, (*THC)[i]->GetPos().x()/mm,
				   (*THC)[i]->GetPos().y()/mm,
				   (*THC)[i]->GetPos().z()/mm);
	  
	    }	  
	  analysis->EndOfEvent(n_hit);	  
	}
      
      GammaRayTelDigitizer * myDM = 
	(GammaRayTelDigitizer*)fDM->FindDigitizerModule( "GammaRayTelDigitizer" );
      myDM->Digitize();
      
#ifdef G4STORE_DATA
      // the whole block is needed only when outfile is active; protect block to avoid
      //     compilations warnings from gcc4.6, Gunter Folger

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
	  
	  (*outFile) << std::setw(7) << event_id << " " << NStrip << 
	   	" " << NPlane << " " << IsX << " " << G4endl;	
	}
      }
#endif      
    }
}













