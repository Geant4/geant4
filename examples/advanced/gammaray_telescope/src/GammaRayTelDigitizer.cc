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
// $Id: GammaRayTelDigitizer.cc,v 1.3 2001-12-04 11:40:28 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTelDigitizer  ------
//           by   F.Longo, R.Giannitrapani & G.Santin (13 nov 2000)
//
// ************************************************************

#include "GammaRayTelDigitizer.hh"
#include "GammaRayTelDigi.hh"
#include "GammaRayTelDigitizerMessenger.hh"

#include "GammaRayTelTrackerHit.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4DigiManager.hh"
#include "G4ios.hh"

//#include "G4CollectionNameVector.hh"
#include "g4std/vector"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelDigitizer::GammaRayTelDigitizer(G4String name)
  :G4VDigitizerModule(name)
{


  G4String colName = "DigitsCollection";
  collectionName.push_back(colName);
  Energy_Threshold = 120.*keV;

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
  	("TrackerDigitizer","DigitsCollection"); // to create the Digi Collection
  
  G4DigiManager* DigiMan = G4DigiManager::GetDMpointer();
  
  G4int THCID; 


  // hits collections IDs
 
  THCID = DigiMan->GetHitsCollectionID("TrackerCollection");

  // Hits collection
  
  GammaRayTelTrackerHitsCollection* THC = 0;

  
  THC = (GammaRayTelTrackerHitsCollection*)
		 (DigiMan->GetHitsCollection(THCID));
  
  
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
	      DigitsCollection->insert(Digi);	
	    }  
	}
      
      G4cout << "Number of tracker digits in this event =  " 
	     << DigitsCollection->entries()  
	// << " " << DigitsCollection->GetName()  
	// << " " << DigitsCollection->GetDMname() 
	     << G4endl;
      
      
      StoreDigiCollection(DigitsCollection);

      G4int DCID = -1;
      if(DCID<0)
	{ 
	  //	  DigiMan->List();
	  DCID = DigiMan->GetDigiCollectionID("TrackerDigitizer/DigitsCollection");
	}
      
      
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....










