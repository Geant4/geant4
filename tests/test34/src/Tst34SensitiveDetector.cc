#include "Tst34SensitiveDetector.hh"
#include "Tst34HitsCollection.hh"
#include "Tst34Hit.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include <iostream>

Tst34SensitiveDetector::Tst34SensitiveDetector(G4String name, Tst34DetectorConstruction* det):
G4VSensitiveDetector(name), G4VGFlashSensitiveDetector(), Detector(det)
{
	//@@@@ xN08SensitiveDetector:: evtl name im constructor des G4VGFlashSensitiveDetector ?
	G4String caloname="Tst34Collection";
	collectionName.insert(caloname);
}

Tst34SensitiveDetector::~Tst34SensitiveDetector() {}

void Tst34SensitiveDetector::Initialize(G4HCofThisEvent*HCE)
{
	cout<<"::Initializing the sensitive detector"<<endl;
	static G4int HCID = -1;
	if(HCID<0){ HCID = GetCollectionID(0); }
	HCE->AddHitsCollection( HCID, caloHitsCollection );
	caloHitsCollection=new 
	Tst34HitsCollection(SensitiveDetectorName,collectionName[0]); // first collection
}

void Tst34SensitiveDetector::EndOfEvent(G4HCofThisEvent*HCE)
{
	if (HCE);	
	
}

G4bool Tst34SensitiveDetector::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
{	//cout<<"This is ProcessHits"<<endl;
	G4double e=aStep->GetTotalEnergyDeposit();
	if(e<=0.)return false;
	
	
	Tst34Hit* caloHit=new Tst34Hit();
	caloHit->SetEdep(e);
	caloHit->SetPos(aStep->GetPreStepPoint()->GetPosition());
	caloHitsCollection->insert(caloHit);
	if (ROhist); 
//	G4VPhysicalVolume* physVol = theTouchable->GetVolume();
	
	return true;
}


G4bool Tst34SensitiveDetector::ProcessHits(G4GFlashSpot*aSpot ,G4TouchableHistory* ROhist)
{	//cout<<"This is ProcessHits"<<endl;
	G4double e=aSpot->GetEnergySpot()->GetEnergy();
	if(e<=0.)return false;
	//@@@@   Tst34SensitiveDetector: vielleicht besser pointer ?
	
	
	Tst34Hit* caloHit=new Tst34Hit();
	caloHit->SetEdep(e);
	caloHit->SetPos(aSpot->GetEnergySpot()->GetPosition());
	caloHitsCollection->insert(caloHit);
	if (ROhist); 
	
	return true;
}
