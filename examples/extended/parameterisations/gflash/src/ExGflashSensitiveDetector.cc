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
// Created by Joanna Weng 26.11.2004
#include "ExGflashSensitiveDetector.hh"
#include "ExGflashHitsCollection.hh"
#include "ExGflashHit.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include <iostream>
using namespace std;
//WARNING :  You have to use also  G4VGFlashSensitiveDetector() as base class
ExGflashSensitiveDetector::ExGflashSensitiveDetector(G4String name, ExGflashDetectorConstruction* det):
G4VSensitiveDetector(name), G4VGFlashSensitiveDetector(), Detector(det)
{
	//@@@@ xN08SensitiveDetector:: evtl name im constructor des G4VGFlashSensitiveDetector ?
	G4String caloname="ExGflashCollection";
	collectionName.insert(caloname);
}

ExGflashSensitiveDetector::~ExGflashSensitiveDetector() {}

void ExGflashSensitiveDetector::Initialize(G4HCofThisEvent*HCE)
{
	G4cout<<"::Initializing the sensitive detector"<<G4endl;
	static G4int HCID = -1;
	if(HCID<0){ HCID = GetCollectionID(0); }
	caloHitsCollection=new 
	ExGflashHitsCollection(SensitiveDetectorName,collectionName[0]); // first collection
	HCE->AddHitsCollection( HCID, caloHitsCollection );
}

void ExGflashSensitiveDetector::EndOfEvent(G4HCofThisEvent*HCE)
{
	if (HCE){;}	
	// No Hadronic
	//    cout<<"IEndOfEvent  " << HCID << endl;
}

G4bool ExGflashSensitiveDetector::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
{	//cout<<"This is ProcessHits"<<endl;
	G4double e=aStep->GetTotalEnergyDeposit();
	if(e<=0.)return false;
	
	G4TouchableHistory* theTouchable
	= (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());

	// enrgy deposited -> make Hit
	//const G4VPhysicalVolume* physVol= aStep->GetPreStepPoint()->GetPhysicalVolume();
	//G4TouchableHistory* theTouchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
	ExGflashHit* caloHit=new ExGflashHit();
	caloHit->SetEdep(e);
	caloHit->SetPos(aStep->GetPreStepPoint()->GetPosition());
	caloHitsCollection->insert(caloHit);
	if (ROhist){;} 
	G4VPhysicalVolume* physVol = theTouchable->GetVolume();
	G4int crystalnum=0;
	for(int i=0;i<100;i++) //@@@@@@@ ExGflashSensitiveDetector:vorsichty
	{
		if(physVol == Detector->GetCristal(i))   crystalnum= i;
	}
	caloHit->SetCrystalNum(crystalnum);
	
	return true;
}

// Separate GFLASH interface
G4bool ExGflashSensitiveDetector::ProcessHits(G4GFlashSpot*aSpot ,G4TouchableHistory* ROhist)
{	//cout<<"This is ProcessHits GFLASH"<<endl;
	G4double e=aSpot->GetEnergySpot()->GetEnergy();
	if(e<=0.)return false;
	
	G4VPhysicalVolume* pCurrentVolume = aSpot->GetTouchableHandle()->GetVolume();
	
	ExGflashHit* caloHit=new ExGflashHit();
	caloHit->SetEdep(e);
	caloHit->SetPos(aSpot->GetEnergySpot()->GetPosition());
	caloHitsCollection->insert(caloHit);
	if (ROhist){;} 
	//cout <<pCurrentVolume->GetName()   << endl;
	G4int crystalnum=0;
	for(int i=0;i<100;i++) //@@@@@@@ ExGflashSensitiveDetector:vorsichty
	{
		if(pCurrentVolume  == Detector->GetCristal(i))   crystalnum= i;
	}
	caloHit->SetCrystalNum(crystalnum);
	
	return true;
}
