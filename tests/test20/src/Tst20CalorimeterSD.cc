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
// $Id: Tst20CalorimeterSD.cc,v 1.5 2007-11-09 18:33:00 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 


#include "Tst20CalorimeterSD.hh"

#include "Tst20CalorHit.hh"
#include "Tst20DetectorConstruction.hh"

#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
  
#include "G4ios.hh"


Tst20CalorimeterSD::Tst20CalorimeterSD(G4String name,	 
				       Tst20DetectorConstruction* det) : G4VSensitiveDetector(name),
									 detector(det)
{
  collectionName.insert("CalCollection");
  hitID = new G4int[500];
}


Tst20CalorimeterSD::~Tst20CalorimeterSD()
{
  delete [] hitID;
}


void Tst20CalorimeterSD::Initialize(G4HCofThisEvent*)
{
  collection = new Tst20CalorHitsCollection(SensitiveDetectorName,collectionName[0]); 
  for (G4int j=0; j<1; j++) {hitID[j] = -1;};
}


G4bool Tst20CalorimeterSD::ProcessHits(G4Step* step,G4TouchableHistory*)
{
  G4double energyDeposit = step->GetTotalEnergyDeposit();

  G4double length = 0.;
  
  if (step->GetTrack()->GetDefinition()->GetPDGCharge() != 0.)
    {
      length = step->GetStepLength();
    }

  if ((energyDeposit == 0.) && (length == 0.)) return false;      

  G4TouchableHistory* theTouchable = (G4TouchableHistory*)(step->GetPreStepPoint()->GetTouchable());
    
  G4VPhysicalVolume* physVol = theTouchable->GetVolume(); 
  //theTouchable->MoveUpHistory();

  G4int number = 0 ;

  if (hitID[number] == -1)
    { 
      Tst20CalorHit* calHit = new Tst20CalorHit();
      if (physVol == detector->GetAbsorber()) calHit->AddEnergyDeposit(energyDeposit,length);
      hitID[number] = collection->insert(calHit) - 1;
      if (verboseLevel > 0)
        G4cout << " New Calorimeter Hit on Tst20: " << number << G4endl;
    }
  else
    { 
      if (physVol == detector->GetAbsorber())
	{
	  (*collection)[hitID[number]]->AddEnergyDeposit(energyDeposit,length);
	}
      if (verboseLevel > 0)
	{
	  G4cout << " Energy added to Tst20: " << number << G4endl; 
	}
    }
    
  return true;
}


void Tst20CalorimeterSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  static G4int HCID = -1;
  if (HCID < 0)
    { 
      HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); 
    }
  HCE->AddHitsCollection(HCID,collection);
}


void Tst20CalorimeterSD::clear()
{} 



void Tst20CalorimeterSD::PrintAll()
{} 


