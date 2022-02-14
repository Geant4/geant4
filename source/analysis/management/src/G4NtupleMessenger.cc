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

using namespace G4Analysis;
using std::to_string;

namespace {

//_____________________________________________________________________________
void WrongParametersWarning(
  const G4String& commandName, std::size_t got, std::size_t expected,
  std::string_view className)
{
  Warn(
    "Got wrong number of \"" + commandName + "\" parameters: " +
    to_string(got) + " instead of " + to_string(expected) + " expected",
    className, "SetNewValue");
}

}

//_____________________________________________________________________________
G4NtupleMessenger::G4NtupleMessenger(G4VAnalysisManager* manager)
  : G4UImessenger(),
    fManager(manager)
{
  fNtupleDir = std::make_unique<G4UIdirectory>("/analysis/ntuple/");
  fNtupleDir->SetGuidance("ntuple control");

  SetActivationCmd();
  SetActivationToAllCmd();
  SetFileNameCmd();
  SetFileNameToAllCmd();
}

//_____________________________________________________________________________
G4NtupleMessenger::~G4NtupleMessenger() = default;

//
// public functions
//

//_____________________________________________________________________________
void G4NtupleMessenger::SetActivationCmd()
{
  auto ntupleId = new G4UIparameter("NtupleId", 'i', false);
  ntupleId->SetGuidance("Ntuple id");
  ntupleId->SetParameterRange("NtupleId>=0");

  auto ntupleActivation = new G4UIparameter("NtupleActivation", 's', true);
  ntupleActivation->SetGuidance("Ntuple activation");
  ntupleActivation->SetDefaultValue("none");

  fSetActivationCmd = std::make_unique<G4UIcommand>("/analysis/ntuple/setActivation", this);
  G4String guidance("Set activation for the ntuple of given id");

  fSetActivationCmd->SetGuidance(guidance);
  fSetActivationCmd->SetParameter(ntupleId);
  fSetActivationCmd->SetParameter(ntupleActivation);
  fSetActivationCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//_____________________________________________________________________________
void G4NtupleMessenger::SetActivationToAllCmd()
{
  fSetActivationAllCmd
    = std::make_unique<G4UIcmdWithABool>("/analysis/ntuple/setActivationToAll", this);
  G4String guidance("Set activation to all ntuples");
  fSetActivationAllCmd->SetGuidance(guidance);
  fSetActivationAllCmd->SetParameterName("AllNtupleActivation",false);
}

//_____________________________________________________________________________
void G4NtupleMessenger::SetFileNameCmd()
{
  auto ntupleId = new G4UIparameter("NtupleId", 'i', false);
  ntupleId->SetGuidance("Ntuple id");
  ntupleId->SetParameterRange("NtupleId>=0");

  auto ntupleFileName = new G4UIparameter("NtupleFileName", 's', true);
  ntupleFileName->SetGuidance("Ntuple file name");
  ntupleFileName->SetDefaultValue("none");

  fSetFileNameCmd = std::make_unique<G4UIcommand>("/analysis/ntuple/setFileName", this);
  G4String guidance("Set file name for the ntuple of given id");

  fSetFileNameCmd->SetGuidance(guidance);
  fSetFileNameCmd->SetParameter(ntupleId);
  fSetFileNameCmd->SetParameter(ntupleFileName);
  fSetFileNameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//_____________________________________________________________________________
void G4NtupleMessenger::SetFileNameToAllCmd()
{
  fSetFileNameAllCmd
    = std::make_unique<G4UIcmdWithAString>("/analysis/ntuple/setFileNameToAll", this);
  G4String guidance("Set file name to all ntuples");
  fSetFileNameAllCmd->SetGuidance(guidance);
  fSetFileNameAllCmd->SetParameterName("AllNtupleFileName",false);
}

//
// public methods
//

//_____________________________________________________________________________
void G4NtupleMessenger::SetNewValue(G4UIcommand* command, G4String newValues)
{
  if ( command == fSetActivationCmd.get() ) {
    // tokenize parameters in a vector
    std::vector<G4String> parameters;
    G4Analysis::Tokenize(newValues, parameters);
    // check consistency
    if ( parameters.size() == command->GetParameterEntries() ) {
      auto counter = 0;
      auto id = G4UIcommand::ConvertToInt(parameters[counter++]);
      auto activation = G4UIcommand::ConvertToBool(parameters[counter++]);
      fManager->SetNtupleActivation(id, activation);
    }
    else {
      // Should never happen but let's check anyway for consistency
      WrongParametersWarning(command->GetCommandName(),
        parameters.size(),command->GetParameterEntries(), fkClass);
    }
  }
  else if ( command == fSetActivationAllCmd.get() ) {
    auto activation = fSetActivationAllCmd->GetNewBoolValue(newValues);
    fManager->SetNtupleActivation(activation);
  }
  else if ( command == fSetFileNameCmd.get() ) {
    // tokenize parameters in a vector
    std::vector<G4String> parameters;
    G4Analysis::Tokenize(newValues, parameters);
    // check consistency
    if ( parameters.size() == command->GetParameterEntries() ) {
      auto counter = 0;
      auto id = G4UIcommand::ConvertToInt(parameters[counter++]);
      auto fileName = parameters[counter++];
      fManager->SetNtupleFileName(id, fileName);
    }
    else {
      // Should never happen but let's check anyway for consistency
      WrongParametersWarning(command->GetCommandName(),
        parameters.size(),command->GetParameterEntries(), fkClass);
    }
  }
  else if ( command == fSetFileNameAllCmd.get() ) {
    auto fileName = newValues;
    fManager->SetNtupleFileName(fileName);
  }
}
