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
#include "Tst34SensitiveDetector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include <iostream>

Tst34SensitiveDetector::Tst34SensitiveDetector(G4String name)
  : G4VSensitiveDetector(name), G4VGFlashSensitiveDetector()
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
                                           G4TouchableHistory*)
{
  G4double e=aStep->GetTotalEnergyDeposit();
  if(e<=0.)return false;

  Tst34Hit* caloHit=new Tst34Hit();
  caloHit->SetEdep(e);
  caloHit->SetPos(aStep->GetPreStepPoint()->GetPosition());
  caloHitsCollection->insert(caloHit);

  return true;
}

G4bool Tst34SensitiveDetector::ProcessHits(G4GFlashSpot*aSpot,
                                           G4TouchableHistory*)
{
  G4double e=aSpot->GetEnergySpot()->GetEnergy();
  if(e<=0.) return false;

  Tst34Hit* caloHit=new Tst34Hit();
  caloHit->SetEdep(e);
  caloHit->SetPos(aSpot->GetEnergySpot()->GetPosition());
  caloHitsCollection->insert(caloHit);

  return true;
}
