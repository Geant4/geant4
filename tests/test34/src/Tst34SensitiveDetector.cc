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
#include "Tst34SensitiveDetector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include <iostream>

Tst34SensitiveDetector::Tst34SensitiveDetector(G4String name,
                                               Tst34DetectorConstruction* det)
  : G4VSensitiveDetector(name), G4VGFlashSensitiveDetector(), Detector(det)
{
  G4String caloname="Tst34Collection";
  collectionName.insert(caloname);
}

Tst34SensitiveDetector::~Tst34SensitiveDetector() {}

void Tst34SensitiveDetector::Initialize(G4HCofThisEvent*)
{
  G4cout << "::Initializing the sensitive detector" << G4endl;
  caloHitsCollection =
    new Tst34HitsCollection(SensitiveDetectorName,collectionName[0]);
}

void Tst34SensitiveDetector::EndOfEvent(G4HCofThisEvent*HCE)
{
  static G4int HCID = -1;
  if(HCID<0){ HCID = GetCollectionID(0); }
  HCE->AddHitsCollection( HCID, caloHitsCollection );
}

G4bool Tst34SensitiveDetector::ProcessHits(G4Step* aStep,
                                           G4TouchableHistory* ROhist)
{
  G4double e=aStep->GetTotalEnergyDeposit();
  if(e<=0.)return false;

  Tst34Hit* caloHit=new Tst34Hit();
  caloHit->SetEdep(e);
  caloHit->SetPos(aStep->GetPreStepPoint()->GetPosition());
  caloHitsCollection->insert(caloHit);
  if (ROhist); 
//  G4VPhysicalVolume* physVol = theTouchable->GetVolume();

  return true;
}

G4bool Tst34SensitiveDetector::ProcessHits(G4GFlashSpot*aSpot,
                                           G4TouchableHistory* ROhist)
{
  G4double e=aSpot->GetEnergySpot()->GetEnergy();
  if(e<=0.) return false;

  Tst34Hit* caloHit=new Tst34Hit();
  caloHit->SetEdep(e);
  caloHit->SetPos(aSpot->GetEnergySpot()->GetPosition());
  caloHitsCollection->insert(caloHit);
  if (ROhist); 

  return true;
}
