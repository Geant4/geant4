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
// $Id$

// Author: Ivana Hrivnacova, 05/05/2015  (ivana@ipno.in2p3.fr)

#include "G4NtupleMessenger.hh"
#include "G4VAnalysisManager.hh"
#include "G4AnalysisUtilities.hh"

#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithABool.hh"

#include <iostream>

using namespace G4Analysis;

//_____________________________________________________________________________
G4NtupleMessenger::G4NtupleMessenger(G4VAnalysisManager* manager)
  : G4UImessenger(),
    fManager(manager),
    fSetActivationCmd(nullptr),
    fSetActivationAllCmd(nullptr)
{
  fNtupleDir = G4Analysis::make_unique<G4UIdirectory>("/analysis/ntuple/");
  fNtupleDir->SetGuidance("ntuple control");
  
  SetActivationCmd();
  SetActivationToAllCmd();
}

//_____________________________________________________________________________
G4NtupleMessenger::~G4NtupleMessenger()
{
}

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

  fSetActivationCmd = G4Analysis::make_unique<G4UIcommand>("/analysis/ntuple/setActivation", this);
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
    = G4Analysis::make_unique<G4UIcmdWithABool>("/analysis/ntuple/setActivationToAll", this);
  G4String guidance("Set activation to all ntuples");
  fSetActivationAllCmd->SetGuidance(guidance);
  fSetActivationAllCmd->SetParameterName("AllNtupleActivation",false);
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
    if ( G4int(parameters.size()) == command->GetParameterEntries() ) {
      auto counter = 0;
      auto id = G4UIcommand::ConvertToInt(parameters[counter++]);
      auto activation = G4UIcommand::ConvertToBool(parameters[counter++]);
      fManager->SetNtupleActivation(id, activation);     
    }
    else {
      // Should never happen but let's check anyway for consistency
      G4ExceptionDescription description;
      description 
        << "Got wrong number of \"" << command->GetCommandName() 
        << "\" parameters: " << parameters.size()
        << " instead of " << command->GetParameterEntries() 
        << " expected" << G4endl;
      G4Exception("G4NtupleMessenger::SetNewValue",
                  "Analysis_W013", JustWarning, description);
    }  
  }
  else if ( command == fSetActivationAllCmd.get() ) {
    auto activation = fSetActivationAllCmd->GetNewBoolValue(newValues);
    fManager->SetNtupleActivation(activation);
  }  
}  
