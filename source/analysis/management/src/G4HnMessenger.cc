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
#include "G4AnalysisMessengerHelper.hh"

#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithABool.hh"

using namespace G4Analysis;

#include <iostream>

G4HnMessenger::G4HnMessenger(G4HnManager& manager)
  : G4UImessenger(),
    fManager(manager),
    fHelper(nullptr),
    fSetHnAsciiCmd(nullptr), 
    fSetHnActivationCmd(nullptr),
    fSetHnActivationAllCmd(nullptr),
    fSetHnPlottingCmd(nullptr),
    fSetHnPlottingAllCmd(nullptr)
{ 
  G4String hnType = fManager.GetHnType();
  hnType.toLower();
  fHelper = G4Analysis::make_unique<G4AnalysisMessengerHelper>(hnType);

  SetHnAsciiCmd();
  SetHnActivationCmd();
  SetHnActivationToAllCmd();
  SetHnPlottingCmd();
  SetHnPlottingToAllCmd();
}

//_____________________________________________________________________________
G4HnMessenger::~G4HnMessenger()
{}

//
// private functions
//

//_____________________________________________________________________________
void G4HnMessenger::SetHnAsciiCmd()
{
  fSetHnAsciiCmd
    = G4Analysis::make_unique<G4UIcmdWithAnInteger>(fHelper->Update("/analysis/HNTYPE_/setAscii"), this);
  fSetHnAsciiCmd->SetGuidance(
    fHelper->Update("Print NDIM_D LOBJECT of given id on ascii file."));

  fSetHnAsciiCmd->SetParameterName("id",false);
  fSetHnAsciiCmd->SetRange("id>=0");
  fSetHnAsciiCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//_____________________________________________________________________________
void G4HnMessenger::SetHnActivationCmd()
{
  auto hnId = new G4UIparameter("id", 'i', false);
  hnId->SetGuidance(fHelper->Update("OBJECT id"));
  hnId->SetParameterRange("id>=0");

  auto hnActivation = new G4UIparameter("hnActivation", 's', true);
  hnActivation->SetGuidance(fHelper->Update("OBJECT activation"));
  hnActivation->SetDefaultValue("none");

  fSetHnActivationCmd 
    = G4Analysis::make_unique<G4UIcommand>(fHelper->Update("/analysis/HNTYPE_/setActivation"), this);
  fSetHnActivationCmd->SetGuidance(
      fHelper->Update("Set activation for the NDIM_D LOBJECT of given id"));
  fSetHnActivationCmd->SetParameter(hnId);
  fSetHnActivationCmd->SetParameter(hnActivation);
  fSetHnActivationCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  

//_____________________________________________________________________________
void G4HnMessenger::SetHnActivationToAllCmd()
{
  fSetHnActivationAllCmd 
    = G4Analysis::make_unique<G4UIcmdWithABool>(fHelper->Update("/analysis/HNTYPE_/setActivationToAll"), this);
  fSetHnActivationAllCmd->SetGuidance(
    fHelper->Update("Set activation to all NDIM_D LOBJECTs"));
  fSetHnActivationAllCmd->SetParameterName("Activation",false);
}  
  
//_____________________________________________________________________________
void G4HnMessenger::SetHnPlottingCmd()
{
  auto hnId = new G4UIparameter("id", 'i', false);
  hnId->SetGuidance(fHelper->Update("OBJECT id"));
  hnId->SetParameterRange("id>=0");

  auto hnPlotting = new G4UIparameter("hnPlotting", 's', true);
  hnPlotting->SetGuidance(fHelper->Update("(In)Activate OBJECT plotting"));
  hnPlotting->SetDefaultValue("none");

  fSetHnPlottingCmd 
    = G4Analysis::make_unique<G4UIcommand>(fHelper->Update("/analysis/HNTYPE_/setPlotting"), this);
  fSetHnPlottingCmd->SetGuidance(
      fHelper->Update("(In)Activate batch plotting of the NDIM_D LOBJECT of given id"));
  fSetHnPlottingCmd->SetParameter(hnId);
  fSetHnPlottingCmd->SetParameter(hnPlotting);
  fSetHnPlottingCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  

//_____________________________________________________________________________
void G4HnMessenger::SetHnPlottingToAllCmd()
{
  fSetHnPlottingAllCmd 
    = G4Analysis::make_unique<G4UIcmdWithABool>(fHelper->Update("/analysis/HNTYPE_/setPlottingToAll"), this);
  fSetHnPlottingAllCmd->SetGuidance(
    fHelper->Update("(In)Activate batch plotting of all NDIM_D LOBJECTs"));
  fSetHnPlottingAllCmd->SetParameterName("Plotting",false);
}  

//
// public methods
//

//_____________________________________________________________________________
void G4HnMessenger::SetNewValue(G4UIcommand* command, G4String newValues)
{
  if ( command == fSetHnAsciiCmd.get() ) { 
    auto id = fSetHnAsciiCmd->GetNewIntValue(newValues);
    fManager.SetAscii(id, true); 
  }      
  else if ( command == fSetHnActivationCmd.get() ) { 
    // tokenize parameters in a vector
    std::vector<G4String> parameters;
    G4Analysis::Tokenize(newValues, parameters);
    // check consistency
    if ( G4int(parameters.size()) == command->GetParameterEntries() ) {
      auto counter = 0;
      auto id = G4UIcommand::ConvertToInt(parameters[counter++]);
      auto activation = G4UIcommand::ConvertToBool(parameters[counter++]);
      fManager.SetActivation(id, activation);     
    }
    else {
      // Should never happen but let's check anyway for consistency
      fHelper->WarnAboutParameters(command, parameters.size());
    }  
  }
  else if ( command == fSetHnActivationAllCmd.get() ) { 
    auto activation = fSetHnActivationAllCmd->GetNewBoolValue(newValues);
    fManager.SetActivation(activation);
  }  
  else if ( command == fSetHnPlottingCmd.get() ) { 
    // tokenize parameters in a vector
    std::vector<G4String> parameters;
    G4Analysis::Tokenize(newValues, parameters);
    // check consistency
    if ( G4int(parameters.size()) == command->GetParameterEntries() ) {
      auto counter = 0;
      auto id = G4UIcommand::ConvertToInt(parameters[counter++]);
      auto activation = G4UIcommand::ConvertToBool(parameters[counter++]);
      fManager.SetPlotting(id, activation);     
    }
    else {
      // Should never happen but let's check anyway for consistency
      fHelper->WarnAboutParameters(command, parameters.size());
    }  
  }
  else if ( command == fSetHnPlottingAllCmd.get() ) { 
    auto activation = fSetHnPlottingAllCmd->GetNewBoolValue(newValues);
    fManager.SetPlotting(activation);
  }  
}  
