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
/// \file ActionInitialization.cc
/// \brief Implementation of the ActionInitialization class

#include "ActionInitialization.hh"

#include "ChemistryWorld.hh"
#include "DetectorConstruction.hh"
#include "InterPulseAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "PulseAction.hh"
#include "RunAction.hh"
#include "StackingAction.hh"
#include "TimeStepAction.hh"

#include "G4DNAChemistryManager.hh"
#include "G4DNAEventScheduler.hh"
#include "G4DNAScavengerMaterial.hh"
#include "G4H2O.hh"
#include "G4Molecule.hh"
#include "G4MoleculeGun.hh"
#include "G4Scheduler.hh"
#include "G4RunManager.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

ActionInitialization::ActionInitialization()
  : G4VUserActionInitialization()
{
  DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ActionInitialization::BuildForMaster() const
{
  SetUserAction(new RunAction());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ActionInitialization::Build() const
{
  PulseAction* pPulseAction = nullptr;
  if (fActivePulse) {  // Le Tuan Anh: define if fActivePulse
    if (fUseInterPulse) {
      pPulseAction =
        new InterPulseAction(fPulseStructure, fUseHistoInput, fPulsePeriod, fNumberOfPulse);
    }
    else {
      pPulseAction = new PulseAction(fPulseStructure, fUseHistoInput);
    }

    pPulseAction->SetPulse(fActivePulse);
    SetUserAction(pPulseAction);
  }

  SetUserAction(new PrimaryGeneratorAction());
  auto pRunAction = new RunAction();
  SetUserAction(pRunAction);

  if (G4DNAChemistryManager::IsActivated()) {
    SetUserAction(new StackingAction());
    const auto* fpDetector = dynamic_cast<const DetectorConstruction*>(
    G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    auto pChemWorld = fpDetector->GetChemistryWorld();
    auto pScavenger = std::make_unique<G4DNAScavengerMaterial>(pChemWorld);
    // To counter Scavenger
    dynamic_cast<G4DNAScavengerMaterial*>(pScavenger.get())->SetCounterAgainstTime();
    G4Scheduler::Instance()->SetScavengerMaterial(std::move(pScavenger));
    auto timeStepAction = new TimeStepAction(pChemWorld, pPulseAction);
    auto eventScheduler = timeStepAction->GetEventScheduler();
    pRunAction->SetEventScheduler(eventScheduler);
    G4Scheduler::Instance()->SetUserAction(timeStepAction);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ActionInitialization::DefineCommands()
{
  // Le Tuan Anh: add new commands:
  fMessenger = std::make_unique<G4GenericMessenger>(this, "/UHDR/pulse/", "Pulse control");
  auto& filenameCmd =
    fMessenger->DeclareMethod("pulseFile", &ActionInitialization::SetPulseStructureInput);
  filenameCmd.SetParameterName("filenamePulse", true);
  filenameCmd.SetDefaultValue("");

  auto& activePulseCmd = fMessenger->DeclareProperty("pulseOn", fActivePulse);
  activePulseCmd.SetParameterName("activatePulse", true);
  activePulseCmd.SetDefaultValue("false");

  auto& filenameCmdHisto =
    fMessenger->DeclareMethod("pulseInHisto", &ActionInitialization::SetPulseStructureHistoInput);
  filenameCmdHisto.SetParameterName("filenameInHisto", true);
  filenameCmdHisto.SetDefaultValue("");

  auto& interPulseCmd = fMessenger->DeclareProperty("multiPulse", fUseInterPulse);
  interPulseCmd.SetParameterName("activateInterPulse", true);
  interPulseCmd.SetDefaultValue("false");

  auto& pulsePeriodCmd =
    fMessenger->DeclareMethodWithUnit("pulsePeriod", "us", &ActionInitialization::SetPulsePeriod);
  pulsePeriodCmd.SetParameterName("pulsePeriod", true);
  pulsePeriodCmd.SetDefaultValue("0");
  pulsePeriodCmd.SetRange("pulsePeriod >= 0");

  auto& nPulseCmd =
    fMessenger->DeclareMethod("numberOfPulse", &ActionInitialization::SetNumberOfPulse);
  nPulseCmd.SetParameterName("numberOfPulse", true);
  nPulseCmd.SetDefaultValue("1");
  nPulseCmd.SetRange("numberOfPulse >= 1");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
void ActionInitialization::SetPulseStructureHistoInput(G4String fn)
{
  // L.T. Anh: set histo input file & activate fUseHistoInput
  fPulseStructure = fn;
  fUseHistoInput = true;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void ActionInitialization::SetPulseStructureInput(G4String fn)
{
  // L.T. Anh: set input file & deactivate fUseHistoInput
  fPulseStructure = fn;
  fUseHistoInput = false;
}