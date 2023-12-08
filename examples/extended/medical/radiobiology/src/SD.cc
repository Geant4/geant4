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
/// \file radiobiology/src/SD.cc
/// \brief Implementation of the RadioBio::SD class

#include "SD.hh"

#include "G4AccumulableManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"

#include "Hit.hh"
#include "Manager.hh"
#include "VRadiobiologicalAccumulable.hh"

namespace RadioBio
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SD::SD(const G4String& name, const G4String& hitsCollectionName) : G4VSensitiveDetector(name)
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fHitsCollection = new RadioBioHitsCollection(SensitiveDetectorName, collectionName[0]);

  // Add this collection in hce
  G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection(hcID, fHitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool SD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  // FIXME why the namespace specification? Should not be necessary...
  RadioBio::Hit* newHit = new RadioBio::Hit();

  // Get the pre-step kinetic energy
  G4double eKinPre = aStep->GetPreStepPoint()->GetKineticEnergy();
  // Get the post-step kinetic energy
  G4double eKinPost = aStep->GetPostStepPoint()->GetKineticEnergy();
  // Get the step average kinetic energy
  G4double eKinMean = (eKinPre + eKinPost) * 0.5;

  const std::vector<const G4Track*>* secondary = aStep->GetSecondaryInCurrentStep();

  size_t SecondarySize = (*secondary).size();
  G4double EnergySecondary = 0.;

  // Get secondary electrons energy deposited
  if (SecondarySize)  // Calculate only secondary particles
  {
    for (size_t numsec = 0; numsec < SecondarySize; numsec++) {
      // Get the PDG code of every secondaty particles in current step
      G4int PDGSecondary = (*secondary)[numsec]->GetDefinition()->GetPDGEncoding();

      if (PDGSecondary == 11)  // Calculate only secondary electrons
      {
        // Calculate the energy deposit of secondary electrons in current step
        EnergySecondary += (*secondary)[numsec]->GetKineticEnergy();
      }
    }
  }

  // Update the hit with the necessary quantities
  newHit->SetTrackID(aStep->GetTrack()->GetTrackID());
  newHit->SetPartType(aStep->GetTrack()->GetParticleDefinition());
  newHit->SetEkinMean(eKinMean);
  newHit->SetDeltaE(aStep->GetTotalEnergyDeposit());
  newHit->SetEinit(aStep->GetPreStepPoint()->GetKineticEnergy());
  newHit->SetTrackLength(aStep->GetStepLength());
  newHit->SetElectronEnergy(EnergySecondary);
  newHit->SetPhysicalVolume(aStep->GetPreStepPoint()->GetPhysicalVolume());
  newHit->SetVoxelIndexes(aStep->GetPreStepPoint()->GetTouchableHandle());

  // Insert the hit in the hitcollection
  fHitsCollection->insert(newHit);

  // Accumulables are accumulated only if calculation is enabled
  for (G4int i = 0; i < G4AccumulableManager::Instance()->GetNofAccumulables(); ++i) {
    // Get the accumulable from the proper manager
    G4VAccumulable* GenAcc = G4AccumulableManager::Instance()->GetAccumulable(i);

    // Get the quantity from the proper manager using the name
    auto q = Manager::GetInstance()->GetQuantity(GenAcc->GetName());

    VRadiobiologicalAccumulable* radioAcc = dynamic_cast<VRadiobiologicalAccumulable*>(GenAcc);

    // If the dynamic_cast did not work, this means that accumulable
    // was not a VRadiobiologicalAccumulable
    if (radioAcc == nullptr) continue;

    // If calculation is enabled, accumulate
    if (q->IsCalculationEnabled()) radioAcc->Accumulate(newHit);
  }

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SD::EndOfEvent(G4HCofThisEvent*)
{
  if (verboseLevel > 1) {
    G4int nofHits = fHitsCollection->entries();
    G4cout << G4endl << "-------->Hits Collection: in this event they are " << nofHits
           << " hits in the detector slices: " << G4endl;
    for (G4int i = 0; i < nofHits; i++)
      (*fHitsCollection)[i]->Print();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace RadioBio