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
// $Id: G4HnMessenger.cc 66310 2012-12-17 11:56:35Z ihrivnac $

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#include "G4HnMessenger.hh"
#include "G4HnManager.hh"
#include "G4AnalysisUtilities.hh"

#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithABool.hh"

#include <iostream>

namespace {

void Exception(G4UIcommand* command, G4int nofParameters)
{
  G4ExceptionDescription description;
  description 
    << "Got wrong number of \"" << command->GetCommandName() 
    << "\" parameters: " << nofParameters
    << " instead of " << command->GetParameterEntries() 
    << " expected" << G4endl;
  G4Exception("G4HnMessenger::SetNewValue",
              "Analysis_W013", JustWarning, description);
}

}                  

//_____________________________________________________________________________
G4HnMessenger::G4HnMessenger(G4HnManager* manager)
  : G4UImessenger(),
    fManager(manager),
    fSetHnAsciiCmd(0), 
    fSetHnActivationCmd(0),
    fSetHnActivationAllCmd(0)
{  
  
  SetHnAsciiCmd();
  SetHnActivationCmd();
  SetHnActivationToAllCmd();
}

//_____________________________________________________________________________
G4HnMessenger::~G4HnMessenger()
{
  delete fSetHnAsciiCmd;  
  delete fSetHnActivationCmd;
  delete fSetHnActivationAllCmd;
}

//
// private functions
//

//_____________________________________________________________________________
G4String  G4HnMessenger::GetCmdDirectoryName() const
{
  G4String hnType = fManager->GetHnType();
  hnType.toLower();
  G4String name("/analysis/");
  name.append(hnType);
  name.append("/");
  return name;
}

//_____________________________________________________________________________
G4String  G4HnMessenger::GetHnDescription() const
{
  G4String hnDescription = fManager->GetHnType();
  hnDescription.erase(0);
  hnDescription.append("D");
  return hnDescription;
}

//_____________________________________________________________________________
void G4HnMessenger::SetHnAsciiCmd()
{
  G4String name = GetCmdDirectoryName();
  name.append("setAscii");
  fSetHnAsciiCmd = new G4UIcmdWithAnInteger(name, this);

  G4String guidance("Print ");
  guidance.append(GetHnDescription());
  guidance.append(" histogram of #Id on ascii file.");

  fSetHnAsciiCmd->SetGuidance(guidance);
  fSetHnAsciiCmd->SetParameterName("Id",false);
  fSetHnAsciiCmd->SetRange("Id>=0");
  fSetHnAsciiCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//_____________________________________________________________________________
void G4HnMessenger::SetHnActivationCmd()
{
  G4UIparameter* hnId = new G4UIparameter("idActivation", 'i', false);
  hnId->SetGuidance("Histogram id");
  hnId->SetParameterRange("idActivation>=0");

  G4UIparameter* hnActivation = new G4UIparameter("hnActivation", 's', true);
  hnActivation->SetGuidance("Histogram activation");
  hnActivation->SetDefaultValue("none");

  G4String name = GetCmdDirectoryName();
  name.append("setActivation");
  fSetHnActivationCmd = new G4UIcommand(name, this);
  
  G4String guidance("Set activation for the ");
  guidance.append(GetHnDescription());
  guidance.append(" histogram of #Id");

  fSetHnActivationCmd->SetGuidance(guidance);
  fSetHnActivationCmd->SetParameter(hnId);
  fSetHnActivationCmd->SetParameter(hnActivation);
  fSetHnActivationCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  

//_____________________________________________________________________________
void G4HnMessenger::SetHnActivationToAllCmd()
{
  G4String name = GetCmdDirectoryName();
  name.append("setActivationToAll");
  fSetHnActivationAllCmd = new G4UIcmdWithABool(name, this);
  
  G4String guidance("Set activation to all ");
  guidance.append(GetHnDescription());
  guidance.append(" histograms.");

  fSetHnActivationAllCmd->SetGuidance(guidance);
  fSetHnActivationAllCmd->SetParameterName("Activation",false);
}  
  
//
// public methods
//

//_____________________________________________________________________________
void G4HnMessenger::SetNewValue(G4UIcommand* command, G4String newValues)
{
  if ( command == fSetHnAsciiCmd ) {
    G4int id = fSetHnAsciiCmd->GetNewIntValue(newValues);
    fManager->SetAscii(id, true); 
  }      
  else if ( command == fSetHnActivationCmd ) {
    // tokenize parameters in a vector
    std::vector<G4String> parameters;
    G4Analysis::Tokenize(newValues, parameters);
    // check consistency
    if ( G4int(parameters.size()) == command->GetParameterEntries() ) {
      G4int counter = 0;
      G4int id = G4UIcommand::ConvertToInt(parameters[counter++]);
      G4bool activation = G4UIcommand::ConvertToBool(parameters[counter++]);
      fManager->SetActivation(id, activation);     
    }
    else {
      // Should never happen but let's check anyway for consistency
      Exception(command, parameters.size());
    }  
  }
  else if ( command == fSetHnActivationAllCmd ) {
    G4bool activation = fSetHnActivationAllCmd->GetNewBoolValue(newValues);
    fManager->SetActivation(activation);
  }  
}  
