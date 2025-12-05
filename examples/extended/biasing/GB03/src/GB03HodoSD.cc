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
/// \file GB03HodoSD.cc
/// \brief Implementation of the GB03HodoSD class

#include "GB03HodoSD.hh"

#include "GB03DetectorConstruction.hh"
#include "GB03Run.hh"

#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GB03HodoSD::GB03HodoSD(G4String name) : G4VSensitiveDetector(name)
{
  collectionName.insert("Leaving");  // Paricles leaving the shield
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GB03HodoSD::Initialize(G4HCofThisEvent* hce)
{
  auto rm = G4RunManager::GetRunManager();
  fRun = static_cast<GB03Run*>(rm->GetNonConstCurrentRun());
  auto det = static_cast<const GB03DetectorConstruction*>(rm->GetUserDetectorConstruction());
  fVerbose = det->GetVerboseLevel();
  // Create hits collection
  fHitsCollection = new GB03HodoHitsCollection(SensitiveDetectorName, collectionName[0]);

  // Add this collection in hce
  auto hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection(hcID, fHitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool GB03HodoSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  if (!step->IsFirstStepInVolume()) return false;

  auto track = step->GetTrack();
  auto preStepPoint = step->GetPreStepPoint();
  auto pname = track->GetParticleDefinition()->GetParticleName();
  auto pid = track->GetParticleDefinition()->GetPDGEncoding();
  auto weight = preStepPoint->GetWeight();
  auto Ekin = preStepPoint->GetKineticEnergy() / MeV;
  auto pos = preStepPoint->GetPosition();

  fRun->CountParticle(pname, weight, Ekin);

  auto hodoHit = new GB03HodoHit();
  hodoHit->Set(pid, Ekin, weight, pos);
  fHitsCollection->insert(hodoHit);

  if (fVerbose > 1) {
    G4cout << std::setw(7) << pid << std::setw(14) << pname
           << ", kinetic energy (MeV) = " << std::setw(12) << Ekin
           << ", position (cm) = " << pos / cm << std::setw(10) << "weight = " << std::setw(10)
           << weight << G4endl;
  }
  return true;
}

void GB03HodoSD::PrintAll()
{
  G4cout << GetFullPathName() << " No of Col " << GetNumberOfCollections() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
