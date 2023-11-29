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
/// \file scavenger/src/PrimaryKiller.cc
/// \brief Implementation of the scavenger::PrimaryKiller class

#include "PrimaryKiller.hh"
#include "G4UnitsTable.hh"
#include "G4Event.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4RunManager.hh"

namespace scavenger
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

PrimaryKiller::PrimaryKiller(const G4String &name, const G4int &depth)
  : G4VPrimitiveScorer(name, depth),
    G4UImessenger(),
    fpELossUI(new G4UIcmdWithADoubleAndUnit("/primaryKiller/eLossMin", this)),
    fpAbortEventIfELossUpperThan(
      new G4UIcmdWithADoubleAndUnit("/primaryKiller/eLossMax", this)),
    fpMinKineticE(new G4UIcmdWithADoubleAndUnit("/primaryKiller/minKineticE", this)),
    fpSizeUI(new G4UIcmdWith3VectorAndUnit("/primaryKiller/setSize", this)) {
  fpSizeUI->SetDefaultUnit("um");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PrimaryKiller::SetNewValue(G4UIcommand *command,
                                G4String newValue) {
  if (command == fpELossUI.get()) {
    fELossRange_Min = fpELossUI->GetNewDoubleValue(newValue);
  } else if (command == fpAbortEventIfELossUpperThan.get()) {
    fELossRange_Max =
      fpAbortEventIfELossUpperThan->GetNewDoubleValue(newValue);
  } else if (command == fpSizeUI.get()) {
    fPhantomSize = fpSizeUI->GetNew3VectorValue(newValue);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4bool PrimaryKiller::ProcessHits(G4Step *aStep, G4TouchableHistory *) {
  const G4Track *track = aStep->GetTrack();
  G4ThreeVector pos = aStep->GetPostStepPoint()->GetPosition();
  if (std::abs(pos.x()) > fPhantomSize.getX() / 2 ||
      std::abs(pos.y()) > fPhantomSize.getY() / 2 ||
      std::abs(pos.z()) > fPhantomSize.getZ() / 2) {
    ((G4Track *) track)->SetTrackStatus(fStopAndKill);
    return false;
  }
  if (track->GetTrackID() != 1) {
    return FALSE;
  }
  G4double kineticE = aStep->GetPostStepPoint()->GetKineticEnergy();
  G4double eLoss = aStep->GetPreStepPoint()->GetKineticEnergy()
                   - kineticE;
  if (eLoss == 0.) {
    return FALSE;
  }
  fELoss += eLoss;
  if (fELoss > fELossRange_Max) {
    G4RunManager::GetRunManager()->AbortEvent();
  }
  if (fELoss >= fELossRange_Min || kineticE <= fKineticE_Min) {
    ((G4Track *) track)->SetTrackStatus(fStopAndKill);
  }
  return TRUE;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PrimaryKiller::Initialize(G4HCofThisEvent * /*HCE*/) {
  fELoss = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

}