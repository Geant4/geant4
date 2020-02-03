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
//      ------------ GammaRayTelDigitizer  ------
//           by   F.Longo, R.Giannitrapani & G.Santin (13 nov 2000)
//
// ************************************************************

#include <vector>

#include "GammaRayTelDigitizer.hh"
#include "GammaRayTelDigi.hh"
#include "GammaRayTelDigitizerMessenger.hh"

#include "GammaRayTelTrackerHit.hh"
#include "GammaRayTelCalorimeterHit.hh"
#include "GammaRayTelAnticoincidenceHit.hh"

#include "G4SystemOfUnits.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4DigiManager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelDigitizer::GammaRayTelDigitizer(G4String name)
  :G4VDigitizerModule(name)
{


  G4String colName = "DigitsCollection";
  collectionName.push_back(colName);
  Energy_Threshold = 120.*keV;
  TotalEnergy = 0.;
  ACDThreshold = 15*keV;

//create a messenger for this class
  digiMessenger = new GammaRayTelDigitizerMessenger(this);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelDigitizer::~GammaRayTelDigitizer()
{
  delete digiMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelDigitizer::Digitize()
{
  DigitsCollection = new GammaRayTelDigitsCollection
  	("GammaRayTelDigitizer","DigitsCollection"); // to create the Digi Collection
  
  G4DigiManager* DigiMan = G4DigiManager::GetDMpointer();
  
  G4int THCID; // TrackerCollection
  G4int CHCID; // CalorimeterCollection
  G4int AHCID; // AnticoincidenceCollection

  TotalEnergy = 0.; // at each event
  
  // TKR Hits collection
  
  THCID = DigiMan->GetHitsCollectionID("TrackerCollection");  
  GammaRayTelTrackerHitsCollection* THC = 0; 
  THC = (GammaRayTelTrackerHitsCollection*)
    (DigiMan->GetHitsCollection(THCID));


  // CAL Hits collection

  CHCID = DigiMan->GetHitsCollectionID("CalorimeterCollection");
  GammaRayTelCalorimeterHitsCollection* CHC = 0;
  CHC = (GammaRayTelCalorimeterHitsCollection*)
		 (DigiMan->GetHitsCollection(CHCID));

  // ACD Hits collection

  AHCID = DigiMan->GetHitsCollectionID("AnticoincidenceCollection");  
  GammaRayTelAnticoincidenceHitsCollection* AHC = 0;
  AHC = (GammaRayTelAnticoincidenceHitsCollection*)
		 (DigiMan->GetHitsCollection(AHCID));

  
  if (THC)
    {
      G4int n_hit = THC->entries();

      for (G4int i=0;i<n_hit;i++)
	{
	  G4double ESil = (*THC)[i]->GetEdepSil();
	  G4int NStrip = (*THC)[i]->GetNStrip();
	  G4int NPlane = (*THC)[i]->GetNSilPlane();
	  G4int IsX = (*THC)[i]->GetPlaneType();
	  
       // digi generation only if energy deposit greater than threshold

	  if (ESil>Energy_Threshold)
	    {
	      GammaRayTelDigi* Digi = new GammaRayTelDigi();
	      Digi->SetPlaneNumber(NPlane);
	      Digi->SetPlaneType(IsX);
	      Digi->SetStripNumber(NStrip);
	      Digi->SetDigiType(0);
	      Digi->SetEnergy(0.);
	      DigitsCollection->insert(Digi);	
	    }  
	}
    }
  if (CHC)
    {
      G4int n_hit = CHC->entries();
      
      for (G4int i=0;i<n_hit;i++)
	{
	  TotalEnergy +=(*CHC)[i]->GetEdepCAL();
	}
      // digi generation only if energy deposit greater than 0.
      
      if (TotalEnergy>0.)
	{
	  GammaRayTelDigi* Digi = new GammaRayTelDigi();
	  Digi->SetPlaneNumber(0);
	  Digi->SetPlaneType(0);
	  Digi->SetStripNumber(0);
	  Digi->SetDigiType(1);
	  Digi->SetEnergy(TotalEnergy);
	  DigitsCollection->insert(Digi);	
	}  

    }

  if (AHC)
    {
      G4int n_hit = AHC->entries();
      
      for (G4int i=0;i<n_hit;i++)
	{
	  G4double energy = (*AHC)[i]->GetEdepACD();
	  G4int type = (*AHC)[i]->GetACDTileNumber();
	  
	  // digi generation only if energy deposit greater than 0.
	  
	  if (energy>ACDThreshold)
	    {
	      GammaRayTelDigi* Digi = new GammaRayTelDigi();
	      Digi->SetPlaneNumber(0);
	      Digi->SetPlaneType(0);
	      Digi->SetStripNumber(type);
	      Digi->SetDigiType(2);
	      Digi->SetEnergy(energy);
	      DigitsCollection->insert(Digi);	
	    }  
	  
	}
    }
  
  if (THC||AHC||CHC){
    G4cout << "Number of digits in this event =  " 
	   << DigitsCollection->entries()  
      // << " " << DigitsCollection->GetName()  
      // << " " << DigitsCollection->GetDMname() 
	   << G4endl;
  }
  
  StoreDigiCollection(DigitsCollection);
  
  G4int DCID = -1;
  if(DCID<0)
    { 
      //	  DigiMan->List();
      DCID = DigiMan->GetDigiCollectionID("GammaRayTelDigitizer/DigitsCollection");
    }
  
  
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....










