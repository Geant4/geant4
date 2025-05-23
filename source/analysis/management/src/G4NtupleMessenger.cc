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

// Author: Ivana Hrivnacova, 05/05/2015  (ivana@ipno.in2p3.fr)

#include "G4NtupleMessenger.hh"
#include "G4VAnalysisManager.hh"
#include "G4AnalysisUtilities.hh"

#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"

using namespace G4Analysis;
using std::to_string;

//_____________________________________________________________________________
G4NtupleMessenger::G4NtupleMessenger(G4VAnalysisManager* manager)
  : fManager(manager)
{
  fNtupleDir = std::make_unique<G4UIdirectory>("/analysis/ntuple/");
  fNtupleDir->SetGuidance("ntuple control");

  CreateCmd();
  CreateColumnCmds();
  FinishCmd();
  DeleteCmd();
  SetActivationCmd();
  SetActivationToAllCmd();
  SetFileNameCmd();
  SetFileNameToAllCmd();
  ListCmd();
}

//_____________________________________________________________________________
G4NtupleMessenger::~G4NtupleMessenger() = default;

//
// private functions
//

//_____________________________________________________________________________
void G4NtupleMessenger::AddIdParameter(G4UIcommand& command)
{
  auto ntupleId = new G4UIparameter("NtupleId", 'i', false);
  ntupleId->SetGuidance("Ntuple id");
  ntupleId->SetParameterRange("NtupleId>=0");

  command.SetParameter(ntupleId);
}

//_____________________________________________________________________________
void G4NtupleMessenger::CreateCmd()
{
  fCreateCmd = CreateCommand<G4UIcommand>("create", "Create ntuple");

  auto ntName = new G4UIparameter("name", 's', false);
  ntName->SetGuidance("Ntuple name");
  fCreateCmd->SetParameter(ntName);

  auto ntTitle = new G4UIparameter("title", 's', false);
  ntTitle->SetGuidance("Ntuple title");
  fCreateCmd->SetParameter(ntTitle);
}

//_____________________________________________________________________________
void G4NtupleMessenger::CreateColumnCmds()
{
  std::vector<char> colTypes = {'I', 'F', 'D', 'S'};

  for (auto colType : colTypes ) {
    G4String name = "createColumn";
    G4String guidance = "Create  ntuple column";
    name.insert(6, 1, colType);
    guidance.insert(7, 1, colType);
    auto cmd = CreateCommand<G4UIcmdWithAString>(std::move(name), std::move(guidance));
    fCreateColumnCmds[colType] = std::move(cmd);
  }
}

//_____________________________________________________________________________
void G4NtupleMessenger::FinishCmd()
{
  fFinishCmd = CreateCommand<G4UIcmdWithoutParameter>(
    "finish", "Finish creating ntuple");
}

//_____________________________________________________________________________
void G4NtupleMessenger::DeleteCmd()
{
  fDeleteCmd = CreateCommand<G4UIcommand>(
    "delete", "Delete ntuple with given id");

  // Add Id parameter
  AddIdParameter(*fDeleteCmd);

  auto parKeepSetting = new G4UIparameter("keepSetting", 'b', true);
  G4String guidance =
    "If set true, activation, file name, etc. options will be kept\n"
    "and applied when a new object with the same id is created.";
  parKeepSetting->SetGuidance(guidance.c_str());
  parKeepSetting->SetDefaultValue("false");
  fDeleteCmd->SetParameter(parKeepSetting);
}

//_____________________________________________________________________________
void G4NtupleMessenger::SetActivationCmd()
{
  fSetActivationCmd = CreateCommand<G4UIcommand>(
    "setActivation", "Set activation for the ntuple with given id");

  AddIdParameter(*fSetActivationCmd);

  auto ntupleActivation = new G4UIparameter("NtupleActivation", 'b', true);
  ntupleActivation->SetGuidance("Ntuple activation");
  ntupleActivation->SetDefaultValue(true);
  fSetActivationCmd->SetParameter(ntupleActivation);
}

//_____________________________________________________________________________
void G4NtupleMessenger::SetActivationToAllCmd()
{
  fSetActivationAllCmd = CreateCommand<G4UIcmdWithABool>(
    "setActivationToAll", "Set activation to all ntuples");
  fSetActivationAllCmd->SetParameterName("AllNtupleActivation",false);
}

//_____________________________________________________________________________
void G4NtupleMessenger::SetFileNameCmd()
{
  fSetFileNameCmd = CreateCommand<G4UIcommand>(
    "setFileName", "Set file name for the ntuple with given id");

  AddIdParameter(*fSetFileNameCmd);

  auto ntupleFileName = new G4UIparameter("NtupleFileName", 's', false);
  ntupleFileName->SetGuidance("Ntuple file name");
  fSetFileNameCmd->SetParameter(ntupleFileName);
}

//_____________________________________________________________________________
void G4NtupleMessenger::SetFileNameToAllCmd()
{
  fSetFileNameAllCmd = CreateCommand<G4UIcmdWithAString>(
    "setFileNameToAll", "Set file name to all ntuples");
  fSetFileNameAllCmd->SetParameterName("AllNtupleFileName",false);
}

//_____________________________________________________________________________
void G4NtupleMessenger::ListCmd()
{
  fListCmd = CreateCommand<G4UIcommand>("list", "List all/active ntuples");
  fListCmd->AvailableForStates(G4State_Idle, G4State_GeomClosed, G4State_EventProc);

  auto parOnlyIfActive = new G4UIparameter("onlyIfActive", 'b', true);
  parOnlyIfActive->SetGuidance("Option whether to list only active objects");
  parOnlyIfActive->SetDefaultValue("true");
  fListCmd->SetParameter(parOnlyIfActive);
}

//
// public methods
//

//_____________________________________________________________________________
void G4NtupleMessenger::SetNewValue(G4UIcommand* command, G4String newValues)
{
  // process "All" commands first

  if ( command == fSetActivationAllCmd.get() ) {
    fManager->SetActivation(fSetActivationAllCmd->GetNewBoolValue(newValues));
    return;
  }

  if ( command == fSetFileNameAllCmd.get() ) {
    fManager->SetFileName(newValues);
    return;
  }

  // Tokenize parameters in a vector
  std::vector<G4String> parameters;
  G4Analysis::Tokenize(newValues, parameters);
  // check consistency
  if ( parameters.size() != command->GetParameterEntries() ) {
    // Should never happen but let's check anyway for consistency
    G4Analysis::Warn(
      "Got wrong number of \"" + command->GetCommandName() +
      "\" parameters: " + std::to_string(parameters.size()) +
      " instead of " + std::to_string(command->GetParameterEntries()) + " expected",
      fkClass, "WarnAboutParameters");
    return;
  }

  auto counter = 0;

  // commands without Id parameter

  if ( command == fCreateCmd.get() ) {
    const auto& name = parameters[counter++];
    const auto& title = parameters[counter++];
    fTmpNtupleId = fManager->CreateNtuple(name, title);
    return;
  }

  for (const auto& [colType, checkCommand] : fCreateColumnCmds) {
    if ( command == checkCommand.get() ) {
      const auto& name = parameters[counter++];
      switch (colType) {
        case 'I':
          fManager->CreateNtupleIColumn(fTmpNtupleId, name);
          return;
        case 'F':
          fManager->CreateNtupleFColumn(fTmpNtupleId, name);
          return;
        case 'D':
          fManager->CreateNtupleDColumn(fTmpNtupleId, name);
          return;
        case 'S':
          fManager->CreateNtupleSColumn(fTmpNtupleId, name);
          return;
        default:
          return;
      }
    }
  }

  if ( command == fFinishCmd.get() ) {
    fManager->FinishNtuple(fTmpNtupleId);
    fTmpNtupleId = G4Analysis::kInvalidId;
    return;
  }

  // commands with Id parameter

  auto id = G4UIcommand::ConvertToInt(parameters[counter++]);

  if ( command == fDeleteCmd.get() ) {
    const auto& keepSetting = G4UIcommand::ConvertToBool(parameters[counter++]);
    fManager->DeleteNtuple(id, keepSetting);
    return;
  }

  if ( command == fSetActivationCmd.get() ) {
    fManager->SetNtupleActivation(id, G4UIcommand::ConvertToBool(parameters[counter++]));
    return;
  }

  if ( command == fSetFileNameCmd.get() ) {
    fManager->SetNtupleFileName(id, parameters[counter++]);
    return;
  }

  if ( command == fListCmd.get() ) {
    const auto& onlyIfActive = G4UIcommand::ConvertToBool(parameters[0]);
    fManager->ListNtuple(onlyIfActive);
    return;
  }
}
