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

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#include "G4HnMessenger.hh"
#include "G4HnManager.hh"

#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"

//_____________________________________________________________________________
G4HnMessenger::G4HnMessenger(G4HnManager& manager)
  : fManager(manager),
    fHnType(manager.GetHnType())
{
  SetHnAsciiCmd();
  SetHnActivationCmd();
  SetHnActivationToAllCmd();
  SetHnPlottingCmd();
  SetHnPlottingToAllCmd();
  SetHnFileNameCmd();
  SetHnFileNameToAllCmd();
}

//_____________________________________________________________________________
G4HnMessenger::~G4HnMessenger() = default;

//
// private functions
//

//_____________________________________________________________________________
G4String G4HnMessenger::GetObjectType() const
{
  return (fHnType[0] == 'h') ?
    fHnType.substr(1,1) + "D histogram" :
    fHnType.substr(1,1) + "D profile";
}

//_____________________________________________________________________________
void G4HnMessenger::AddIdParameter(G4UIcommand& command)
{
  auto htId = new G4UIparameter("id", 'i', false);
  htId->SetGuidance("Histogram id");
  htId->SetParameterRange("id>=0");
  command.SetParameter(htId);
}

//_____________________________________________________________________________
void G4HnMessenger::AddOptionParameter(G4UIcommand& command, G4String optionName)
{
  auto param = new G4UIparameter(optionName, 'b', true);
  auto guidance = GetObjectType() + " " + optionName + " option";
  param->SetGuidance(guidance.c_str());
  param->SetDefaultValue("true");
  command.SetParameter(param);
}

//_____________________________________________________________________________
void G4HnMessenger::SetHnAsciiCmd()
{
  fSetAsciiCmd =
    CreateCommand<G4UIcommand>("setAscii", "Print  on ascii file the ");

  AddIdParameter(*fSetAsciiCmd);
  AddOptionParameter(*fSetAsciiCmd, "hnAscii");

}

//_____________________________________________________________________________
void G4HnMessenger::SetHnActivationCmd()
{
  fSetActivationCmd =
    CreateCommand<G4UIcommand>("setActivation", "Set activation to the ");

  AddIdParameter(*fSetActivationCmd);
  AddOptionParameter(*fSetActivationCmd, "hnActivation");
}

//_____________________________________________________________________________
void G4HnMessenger::SetHnActivationToAllCmd()
{
  fSetActivationAllCmd =
    CreateCommand<G4UIcmdWithABool>(
      "setActivationToAll", "Set activation to all");
  fSetActivationAllCmd->SetParameterName("Activation", false);
}

//_____________________________________________________________________________
void G4HnMessenger::SetHnPlottingCmd()
{
  fSetPlottingCmd =
    CreateCommand<G4UIcommand>("setPlotting", "(In)Activate batch plotting of the  ");

  AddIdParameter(*fSetPlottingCmd);
  AddOptionParameter(*fSetPlottingCmd, "hnPlotting");
}

//_____________________________________________________________________________
void G4HnMessenger::SetHnPlottingToAllCmd()
{
  fSetPlottingAllCmd =
    CreateCommand<G4UIcmdWithABool>(
      "setPlottingToAll", "(In)Activate batch plotting of all ");
  fSetPlottingAllCmd->SetParameterName("Plotting", false);
}

//_____________________________________________________________________________
void G4HnMessenger::SetHnFileNameCmd()
{
  fSetFileNameCmd =
    CreateCommand<G4UIcommand>("setFileName", "Set the output file name for the ");

  AddIdParameter(*fSetFileNameCmd);

  auto param = new G4UIparameter("hnFileName", 's', false);
  auto guidance = GetObjectType() + " output file name";
  param->SetGuidance(guidance.c_str());
  fSetFileNameCmd->SetParameter(param);
}

//_____________________________________________________________________________
void G4HnMessenger::SetHnFileNameToAllCmd()
{
  fSetFileNameAllCmd =
    CreateCommand<G4UIcmdWithAString>(
      "setFileNameToAll", "Set output file name for all  ");
   fSetFileNameAllCmd->SetParameterName("FileName", false);
}

//
// public methods
//

//_____________________________________________________________________________
void G4HnMessenger::SetNewValue(G4UIcommand* command, G4String newValues)
{
  // process "All" commands first
  if ( command == fSetActivationAllCmd.get() ) {
    fManager.SetActivation(fSetActivationAllCmd->GetNewBoolValue(newValues));
    return;
  }

  if ( command == fSetPlottingAllCmd.get() ) {
    fManager.SetPlotting(fSetPlottingAllCmd->GetNewBoolValue(newValues));
    return;
  }

  if ( command == fSetFileNameAllCmd.get() ) {
    fManager.SetFileName(newValues);
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
  auto id = G4UIcommand::ConvertToInt(parameters[counter++]);

  if ( command == fSetAsciiCmd.get() ) {
    fManager.SetAscii(id, G4UIcommand::ConvertToBool(parameters[counter++]));
    return;
  }

  if ( command == fSetActivationCmd.get() ) {
    fManager.SetActivation(id, G4UIcommand::ConvertToBool(parameters[counter++]));
    return;
  }

  if ( command == fSetPlottingCmd.get() ) {
    fManager.SetPlotting(id, G4UIcommand::ConvertToBool(parameters[counter++]));
    return;
  }

  if ( command == fSetFileNameCmd.get() ) {
    fManager.SetFileName(id, parameters[counter++]);
    return;
  }
}
