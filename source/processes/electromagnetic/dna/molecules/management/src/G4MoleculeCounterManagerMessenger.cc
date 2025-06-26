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
// Author: Christian Velten (2025)

#include "G4MoleculeCounterManagerMessenger.hh"

#include "G4MoleculeCounterManager.hh"

G4MoleculeCounterManagerMessenger::G4MoleculeCounterManagerMessenger(
  G4MoleculeCounterManager* manager)
  : G4UImessenger(), fpManager(manager)
{
  fpManagerDir = std::make_unique<G4UIdirectory>("/chem/moleculecounters/", false);
  fpManagerDir->SetGuidance("Molecule Counter Manager Commands");
  fpManagerDir->AvailableForStates(G4State_PreInit);

  InitializeCommands();
}

void G4MoleculeCounterManagerMessenger::InitializeCommands()
{
  fpActiveCmd = std::make_unique<G4UIcmdWithABool>("/chem/moleculecounters/active", this);
  fpActiveCmd->AvailableForStates(G4State_PreInit);

  fpResetBeforeEventCmd =
    std::make_unique<G4UIcmdWithABool>("/chem/moleculecounters/resetBeforeEvent", this);
  fpResetBeforeEventCmd->AvailableForStates(G4State_PreInit);

  fpResetBeforeRunCmd =
    std::make_unique<G4UIcmdWithABool>("/chem/moleculecounters/resetBeforeRun", this);
  fpResetBeforeRunCmd->AvailableForStates(G4State_PreInit);

  fpAccumulateIntoMasterCmd =
    std::make_unique<G4UIcmdWithABool>("/chem/moleculecounters/accumulateIntoMaster", this);
  fpAccumulateIntoMasterCmd->AvailableForStates(G4State_PreInit);

  fpVerboseCmd = std::make_unique<G4UIcmdWithAnInteger>("/chem/moleculecounters/verbose", this);
  fpVerboseCmd->SetDefaultValue(0);
  fpVerboseCmd->AvailableForStates(G4State_Idle, G4State_Init, G4State_PreInit, G4State_EventProc,
                                   G4State_GeomClosed);
}

void G4MoleculeCounterManagerMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fpActiveCmd.get()) {
    auto value = fpActiveCmd->ConvertToBool(newValue);
    fpManager->SetIsActive(value);
  }
  else if (command == fpResetBeforeEventCmd.get()) {
    auto value = fpResetBeforeEventCmd->ConvertToBool(newValue);
    fpManager->SetResetCountersBeforeEvent(value);
  }
  else if (command == fpResetBeforeRunCmd.get()) {
    auto value = fpResetBeforeRunCmd->ConvertToBool(newValue);
    fpManager->SetResetCountersBeforeRun(value);
  }
  else if (command == fpAccumulateIntoMasterCmd.get()) {
    auto value = fpAccumulateIntoMasterCmd->ConvertToBool(newValue);
    fpManager->SetAccumulateCounterIntoMaster(value);
  }
  else if (command == fpVerboseCmd.get()) {
    auto value = fpVerboseCmd->ConvertToInt(newValue);
    fpManager->SetVerbosity(value);
  }
}
