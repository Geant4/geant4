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
#include "SiliconPixelSD.hh"

#include "G4String.hh"

SiliconPixelSD::SiliconPixelSD(G4String name)
    : G4VSensitiveDetector("SiliconPixelHitCollection") {
  G4cout << "creating a sensitive detector with name: " << name << G4endl;
  collectionName.insert("SiliconPixelHitCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SiliconPixelSD::~SiliconPixelSD() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SiliconPixelSD::Initialize(G4HCofThisEvent *HCE) {
  fHitCollection = new SiliconPixelHitCollection(GetName(), collectionName[0]);

  if (fHCID < 0)
    fHCID = GetCollectionID(0);
  HCE->AddHitsCollection(fHCID, fHitCollection);

  fTmpHits.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SiliconPixelSD::EndOfEvent(G4HCofThisEvent *) {
  for (auto it = fTmpHits.begin(); it != fTmpHits.end(); ++it)
    fHitCollection->insert(it->second);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool SiliconPixelSD::ProcessHits(G4Step *step, G4TouchableHistory *) {
  G4TouchableHandle touchable = step->GetPreStepPoint()->GetTouchableHandle();

  G4int copyNumCell = touchable->GetVolume(0)->GetCopyNo();
  G4int copyNumSensor = touchable->GetVolume(1)->GetCopyNo();
  int tmp_ID = 1000 * copyNumSensor + copyNumCell;
  if (fTmpHits.find(tmp_ID) == fTmpHits.end()) { // make new hit
    G4String vol_name = touchable->GetVolume(0)->GetName();
    fTmpHits[tmp_ID] =
        new SiliconPixelHit(vol_name, copyNumSensor, copyNumCell);
    G4double hitX = (touchable->GetVolume(1)->GetTranslation().x() +
                     touchable->GetVolume(0)->GetTranslation().x()) /
                    CLHEP::cm;
    G4double hitY = (touchable->GetVolume(1)->GetTranslation().y() +
                     touchable->GetVolume(0)->GetTranslation().y()) /
                    CLHEP::cm;
    G4double hitZ = touchable->GetVolume(1)->GetTranslation().z() / CLHEP::cm;
    fTmpHits[tmp_ID]->SetPosition(hitX, hitY, hitZ); // in cm
  }

  G4double edep = step->GetTotalEnergyDeposit() / CLHEP::keV; // in keV
  G4double edepNonIonizing = step->GetNonIonizingEnergyDeposit() / CLHEP::keV;

  G4double timedep = step->GetPostStepPoint()->GetGlobalTime() / CLHEP::ns;

  fTmpHits[tmp_ID]->AddEdep(edep, timedep);
  fTmpHits[tmp_ID]->AddEdepNonIonizing(edepNonIonizing, timedep);

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
