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
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"
#include "G4DNAChemistryManager.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction() {
  fpDetector = dynamic_cast<const DetectorConstruction *>(
    G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  fRNP = fpDetector->GetNPRadius() / CLHEP::nm;
  fRAbs = fpDetector->GetAbsRadius() / CLHEP::nm;
  fTrackCut = fpDetector->GetTrackingCut() / CLHEP::eV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step *aStep) {
  // Refresh detector parameters (kept to preserve original behavior)
  fpDetector = dynamic_cast<const DetectorConstruction *>(
    G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  fRNP = fpDetector->GetNPRadius() / CLHEP::nm;
  fRAbs = fpDetector->GetAbsRadius() / CLHEP::nm;
  fTrackCut = fpDetector->GetTrackingCut() / CLHEP::eV;

  const auto pre = aStep->GetPreStepPoint();
  const auto post = aStep->GetPostStepPoint();

  const G4ThreeVector &pos = pre->GetPosition();
  const G4ThreeVector &postpos = post->GetPosition();

  const G4double R = pos.mag() / CLHEP::nm;

  G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();

  // Radial energy deposition
  if (fRNP < R) {
    const G4double energy = aStep->GetTotalEnergyDeposit() / CLHEP::joule;
    analysisManager->FillH1(1, R, energy);

    // Angular distribution in the near-equatorial plane
    if (std::abs(pos.z()) < 10 * CLHEP::nm) {
      const G4double theta = std::atan2(pos.y(), pos.x()) / CLHEP::deg;
      if (0 <= theta) {
        analysisManager->FillH2(0, theta, R, energy);
      } else {
        analysisManager->FillH2(0, theta + 360, R, energy);
      }
    }
  }

  // Process tracks crossing NanoParticle boundary
  if (pre->GetPhysicalVolume()->GetName() != "NanoParticle") {
    return;
  }
  if (post->GetStepStatus() != fGeomBoundary) {
    return;
  }

  //*** WARNING: this will kill all incident electrons at the NP surface (as in original code) ***
  G4Track *track = aStep->GetTrack();
  const G4double trackE = track->GetKineticEnergy() / CLHEP::eV;

  if (track->GetTrackID() == 1) {
    if (pos.x() < 0) {
      analysisManager->FillH1(8, trackE);
    } else {
      analysisManager->FillH1(9, trackE);
    }
    if (!G4DNAChemistryManager::IsActivated()) {
      track->SetTrackStatus(fStopAndKill);
    }
    return;
  }

  const G4ThreeVector &dir = post->GetMomentumDirection();
  if (dir.dot(postpos) < 0.0) {
    return;
  }

  if (trackE < fTrackCut) {
    track->SetTrackStatus(fStopAndKill);
    return;
  }

  if (track->GetDefinition()->GetPDGCharge() != 0) {
    analysisManager->FillH1(4, trackE);
  } else {
    analysisManager->FillH1(5, trackE);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......