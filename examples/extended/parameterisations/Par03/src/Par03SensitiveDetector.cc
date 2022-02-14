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
#include "Par03SensitiveDetector.hh"
#include "Par03Hit.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"

Par03SensitiveDetector::Par03SensitiveDetector(G4String aName)
  : G4VSensitiveDetector(aName)
{
  collectionName.insert("hits");
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par03SensitiveDetector::Par03SensitiveDetector(G4String aName, G4int aNumLayers,
                                               G4int aNumRho, G4int aNumPhi)
  : G4VSensitiveDetector(aName)
  , fCellNoZ(aNumLayers)
  , fCellNoRho(aNumRho)
  , fCellNoPhi(aNumPhi)
{
  collectionName.insert("hits");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par03SensitiveDetector::~Par03SensitiveDetector() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par03SensitiveDetector::Initialize(G4HCofThisEvent* aHCE)
{
  fHitsCollection =
    new Par03HitsCollection(SensitiveDetectorName, collectionName[0]);
  if(fHitCollectionID < 0)
  {
    fHitCollectionID =
      G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
  }
  aHCE->AddHitsCollection(fHitCollectionID, fHitsCollection);

  // fill calorimeter hits with zero energy deposition
  for(G4int iphi = 0; iphi < fCellNoPhi; iphi++)
    for(G4int irho = 0; irho < fCellNoRho; irho++)
      for(G4int iz = 0; iz < fCellNoZ; iz++)
      {
        Par03Hit* hit = new Par03Hit();
        fHitsCollection->insert(hit);
      }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool Par03SensitiveDetector::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep == 0.)
    return true;

  G4TouchableHistory* aTouchable =
    (G4TouchableHistory*) (aStep->GetPreStepPoint()->GetTouchable());

  auto hit = RetrieveAndSetupHit(aTouchable);

  // Add energy deposit from G4Step
  hit->AddEdep(edep);

  // Fill time information from G4Step
  // If it's already filled, choose hit with earliest global time
  if(hit->GetTime() == -1 ||
     hit->GetTime() > aStep->GetTrack()->GetGlobalTime())
    hit->SetTime(aStep->GetTrack()->GetGlobalTime());

  // Set hit type to full simulation (only if hit is not already marked as fast
  // sim)
  if(hit->GetType() != 1)
    hit->SetType(0);

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool Par03SensitiveDetector::ProcessHits(const G4FastHit* aHit,
                                           const G4FastTrack* aTrack,
                                           G4TouchableHistory* aTouchable)
{
  G4double edep = aHit->GetEnergy();
  if(edep == 0.)
    return true;

  auto hit = RetrieveAndSetupHit(aTouchable);

  // Add energy deposit from G4FastHit
  hit->AddEdep(edep);

  // Fill time information from G4FastTrack
  // If it's already filled, choose hit with earliest global time
  if(hit->GetTime() == -1 ||
     hit->GetTime() > aTrack->GetPrimaryTrack()->GetGlobalTime())
  {
    hit->SetTime(aTrack->GetPrimaryTrack()->GetGlobalTime());
  }

  // Set hit type to fast simulation (even if hit was already marked as full
  // sim, overwrite it)
  hit->SetType(1);

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par03Hit* Par03SensitiveDetector::RetrieveAndSetupHit(
  G4TouchableHistory* aTouchable)
{
  G4int rhoNo = aTouchable->GetCopyNumber(0);  // cell
  G4int phiNo = aTouchable->GetCopyNumber(1);  // segment
  G4int zNo   = aTouchable->GetCopyNumber(2);  // layer

  std::size_t hitID = fCellNoRho * fCellNoZ * phiNo + fCellNoZ * rhoNo + zNo;

  if(hitID >= fHitsCollection->entries())
  {
    G4Exception(
      "Par03SensitiveDetector::RetrieveAndSetupHit()", "InvalidSetup",
      FatalException,
      "Size of hit collection in Par03SensitiveDetector is smaller than the "
      "number of cells created in Par03DetectorConstruction!");
  }
  Par03Hit* hit = (*fHitsCollection)[hitID];

  if(hit->GetRhoId() < 0)
  {
    hit->SetRhoId(rhoNo);
    hit->SetPhiId(phiNo);
    hit->SetZid(zNo);
    hit->SetLogV(aTouchable->GetVolume(0)->GetLogicalVolume());
    G4AffineTransform transform = aTouchable->GetHistory()->GetTopTransform();
    hit->SetRot(transform.NetRotation());
    transform.Invert();
    hit->SetPos(transform.NetTranslation());
  }
  return hit;
}