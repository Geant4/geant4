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
/// \file B5/src/HadCalorimeterSD.cc
/// \brief Implementation of the B5::HadCalorimeterSD class

#include "HadCalorimeterSD.hh"
#include "HadCalorimeterHit.hh"
#include "Constants.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

namespace B5
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HadCalorimeterSD::HadCalorimeterSD(G4String name)
: G4VSensitiveDetector(name)
{
  collectionName.insert("HadCalorimeterColl");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HadCalorimeterSD::Initialize(G4HCofThisEvent* hce)
{
  fHitsCollection
    = new HadCalorimeterHitsCollection(SensitiveDetectorName,collectionName[0]);
  if (fHCID<0) {
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
  }
  hce->AddHitsCollection(fHCID,fHitsCollection);

  // fill calorimeter hits with zero energy deposition
  for (auto column=0;column<kNofHadColumns;column++) {
    for (auto row=0;row<kNofHadRows;row++) {
      fHitsCollection->insert(new HadCalorimeterHit());
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool HadCalorimeterSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  auto edep = step->GetTotalEnergyDeposit();
  if (edep==0.) return true;

  auto touchable = step->GetPreStepPoint()->GetTouchable();
  auto rowNo = touchable->GetCopyNumber(2);
  auto columnNo = touchable->GetCopyNumber(3);
  auto hitID = kNofHadRows*columnNo+rowNo;
  auto hit = (*fHitsCollection)[hitID];

  // check if it is first touch
  if (hit->GetColumnID()<0) {
    hit->SetColumnID(columnNo);
    hit->SetRowID(rowNo);
    auto depth = touchable->GetHistory()->GetDepth();
    auto transform = touchable->GetHistory()->GetTransform(depth-2);
    transform.Invert();
    hit->SetRot(transform.NetRotation());
    hit->SetPos(transform.NetTranslation());
  }
  // add energy deposition
  hit->AddEdep(edep);

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
