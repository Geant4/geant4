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

// Author: Ivana Hrivnacova, 24/07/2014  (ivana@ipno.in2p3.fr)

#include "G4P1Messenger.hh"
#include "G4VAnalysisManager.hh"
#include "G4AnalysisUtilities.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"

#include <iostream>

using namespace G4Analysis;

//_____________________________________________________________________________
G4P1Messenger::G4P1Messenger(G4VAnalysisManager* manager)
  : G4UImessenger(),
    fManager(manager),
    fHelper(nullptr),
    fDirectory(nullptr),
    fCreateP1Cmd(nullptr),
    fSetP1Cmd(nullptr),
    fSetP1TitleCmd(nullptr), 
    fSetP1XAxisCmd(nullptr), 
    fSetP1YAxisCmd(nullptr),
    fXData()
{  
  fHelper = G4Analysis::make_unique<G4AnalysisMessengerHelper>("p1");

  fDirectory = fHelper->CreateHnDirectory();

  CreateP1Cmd();

  SetP1Cmd();
  fSetP1XCmd = fHelper->CreateSetBinsCommand("x", this);
  fSetP1YCmd = fHelper->CreateSetValuesCommand("y", this);
   
  fSetP1TitleCmd = fHelper->CreateSetTitleCommand(this);
  fSetP1XAxisCmd = fHelper->CreateSetAxisCommand("x", this);
  fSetP1YAxisCmd = fHelper->CreateSetAxisCommand("y", this);
}

//_____________________________________________________________________________
G4P1Messenger::~G4P1Messenger()
{}

//
// private functions
//

//_____________________________________________________________________________
void G4P1Messenger::CreateP1Cmd()
{
  auto p1Name = new G4UIparameter("name", 's', false);
  p1Name->SetGuidance("Profile name (label)");
  
  auto p1Title = new G4UIparameter("title", 's', false);
  p1Title->SetGuidance("Profile title");

  auto p1xNbins0 = new G4UIparameter("xnbins0", 'i', true);
  p1xNbins0->SetGuidance("Number of x-bins (default = 100)");
  p1xNbins0->SetGuidance("Can be reset with /analysis/p1/set command");
  p1xNbins0->SetDefaultValue(100);
  
  auto p1xValMin0 = new G4UIparameter("xvalMin0", 'd', true);
  p1xValMin0->SetGuidance("Minimum x-value, expressed in unit (default = 0.)");
  p1xValMin0->SetGuidance("Can be reset with /analysis/p1/set command");
  p1xValMin0->SetDefaultValue(0.);
  
  auto p1xValMax0 = new G4UIparameter("xvalMax0", 'd', true);
  p1xValMax0->SetGuidance("Maximum x-value, expressed in unit (default = 1.)");
  p1xValMax0->SetGuidance("Can be reset with /analysis/p1/set command");
  p1xValMax0->SetDefaultValue(1.);

  auto p1xValUnit0 = new G4UIparameter("xvalUnit0", 's', true);
  p1xValUnit0->SetGuidance("The unit applied to filled x-values and xvalMin0, xvalMax0");
  p1xValUnit0->SetDefaultValue("none");
  
  auto p1xValFcn0 = new G4UIparameter("xvalFcn0", 's', true);
  G4String fcnxGuidance = "The function applied to filled x-values (log, log10, exp, none).\n";
  fcnxGuidance += "Note that the unit parameter cannot be omitted in this case,\n";
  fcnxGuidance += "but none value should be used insted.";
  p1xValFcn0->SetGuidance(fcnxGuidance);
  p1xValFcn0->SetParameterCandidates("log log10 exp none");
  p1xValFcn0->SetDefaultValue("none");
    
  auto p1xValBinScheme0 = new G4UIparameter("xvalBinScheme0", 's', true);
  G4String binSchemeGuidance = "The binning scheme (linear, log).\n";
  p1xValBinScheme0->SetParameterCandidates("linear log");
  binSchemeGuidance 
    += "Note that the unit and fcn parameters cannot be omitted in this case,\n";
  binSchemeGuidance += "but none value should be used insted.";
  p1xValBinScheme0->SetGuidance(binSchemeGuidance);
  p1xValBinScheme0->SetDefaultValue("linear");
  
  auto p1yValMin0 = new G4UIparameter("yvalMin0", 'd', true);
  p1yValMin0->SetGuidance("Minimum y-value, expressed in unit (default = 0.)");
  p1yValMin0->SetGuidance("Can be reset with /analysis/p1/set command");
  p1yValMin0->SetDefaultValue(0.);
  
  auto p1yValMax0 = new G4UIparameter("yvalMax0", 'd', true);
  p1yValMax0->SetGuidance("Maximum y-value, expressed in unit (default = 1.)");
  p1yValMax0->SetGuidance("Can be reset with /analysis/p1/set command");
  p1yValMax0->SetDefaultValue(1.);

  auto p1yValUnit0 = new G4UIparameter("yvalUnit0", 's', true);
  p1yValUnit0->SetGuidance("The unit applied to filled y-values and yvalMin0, yvalMax0");
  p1yValUnit0->SetDefaultValue("none");
  
  auto p1yValFcn0 = new G4UIparameter("yvalFcn0", 's', true);
  G4String fcnyGuidance = "The function applied to filled y-values (log, log10, exp, none).\n";
  fcnyGuidance += "Note that the unit parameter cannot be omitted in this case,\n";
  fcnyGuidance += "but none value should be used insted.";
  p1yValFcn0->SetGuidance(fcnyGuidance);
  p1yValFcn0->SetParameterCandidates("log log10 exp none");
  p1yValFcn0->SetDefaultValue("none");
  
  fCreateP1Cmd = G4Analysis::make_unique<G4UIcommand>("/analysis/p1/create", this);
  fCreateP1Cmd->SetGuidance("Create 1D profile");
  fCreateP1Cmd->SetParameter(p1Name);
  fCreateP1Cmd->SetParameter(p1Title);
  fCreateP1Cmd->SetParameter(p1xNbins0);
  fCreateP1Cmd->SetParameter(p1xValMin0);
  fCreateP1Cmd->SetParameter(p1xValMax0);
  fCreateP1Cmd->SetParameter(p1xValUnit0);
  fCreateP1Cmd->SetParameter(p1xValFcn0);
  fCreateP1Cmd->SetParameter(p1xValBinScheme0);
  fCreateP1Cmd->SetParameter(p1yValMin0);
  fCreateP1Cmd->SetParameter(p1yValMax0);
  fCreateP1Cmd->SetParameter(p1yValUnit0);
  fCreateP1Cmd->SetParameter(p1yValFcn0);
  fCreateP1Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  


//_____________________________________________________________________________
void G4P1Messenger::SetP1Cmd()
{
  auto p1Id = new G4UIparameter("id", 'i', false);
  p1Id->SetGuidance("Profile id");
  p1Id->SetParameterRange("id>=0");
  
  auto p1xNbins = new G4UIparameter("xnbins", 'i', false);
  p1xNbins->SetGuidance("Number of x-bins");
  
  auto p1xValMin = new G4UIparameter("xvalMin", 'd', false);
  p1xValMin->SetGuidance("Minimum x-value, expressed in unit");
  
  auto p1xValMax = new G4UIparameter("xvalMax", 'd', false);
  p1xValMax->SetGuidance("Maximum x-value, expressed in unit");
  
  auto p1xValUnit = new G4UIparameter("xvalUnit", 's', true);
  p1xValUnit->SetGuidance("The unit applied to filled x-values and xvalMin0, xvalMax0");
  p1xValUnit->SetDefaultValue("none");
 
  auto p1xValFcn = new G4UIparameter("xvalFcn", 's', true);
  p1xValFcn->SetParameterCandidates("log log10 exp none");
  G4String fcnxGuidance = "The function applied to filled x-values (log, log10, exp, none).\n";
  fcnxGuidance += "Note that the unit parameter cannot be omitted in this case,\n";
  fcnxGuidance += "but none value should be used insted.";
  p1xValFcn->SetGuidance(fcnxGuidance);
  p1xValFcn->SetDefaultValue("none");
 
  auto p1xValBinScheme = new G4UIparameter("xvalBinScheme", 's', true);
  G4String binSchemeGuidance = "The binning scheme (linear, log).\n";
  p1xValBinScheme->SetParameterCandidates("linear log");
  binSchemeGuidance 
    += "Note that the unit and fcn parameters cannot be omitted in this case,\n";
  binSchemeGuidance += "but none value should be used insted.";
  p1xValBinScheme->SetGuidance(binSchemeGuidance);
  p1xValBinScheme->SetDefaultValue("linear");
  
  auto p1yValMin = new G4UIparameter("yvalMin", 'd', true);
  p1yValMin->SetGuidance("Minimum y-value, expressed in unit");
  
  auto p1yValMax = new G4UIparameter("yvalMax", 'd', true);
  p1yValMax->SetGuidance("Maximum y-value, expressed in unit");
  
  auto p1yValUnit = new G4UIparameter("yvalUnit", 's', true);
  p1yValUnit->SetGuidance("The unit applied to filled y-values and yvalMin0, yvalMax0");
  p1yValUnit->SetDefaultValue("none");
 
  auto p1yValFcn = new G4UIparameter("yvalFcn", 's', true);
  p1yValFcn->SetParameterCandidates("log log10 exp none");
  G4String fcnyGuidance = "The function applied to filled y-values (log, log10, exp, none).\n";
  fcnyGuidance += "Note that the unit parameter cannot be omitted in this case,\n";
  fcnyGuidance += "but none value should be used insted.";
  p1yValFcn->SetGuidance(fcnyGuidance);
  p1yValFcn->SetDefaultValue("none");
 
  fSetP1Cmd = G4Analysis::make_unique<G4UIcommand>("/analysis/p1/set", this);
  fSetP1Cmd->SetGuidance("Set parameters for the 1D profile of given id:");
  fSetP1Cmd->SetGuidance("  nbins; xvalMin; xvalMax; xunit; xfunction; xbinScheme");
  fSetP1Cmd->SetGuidance("  yvalMin; yvalMax; yunit; yfunction");
  fSetP1Cmd->SetParameter(p1Id);
  fSetP1Cmd->SetParameter(p1xNbins);
  fSetP1Cmd->SetParameter(p1xValMin);
  fSetP1Cmd->SetParameter(p1xValMax);
  fSetP1Cmd->SetParameter(p1xValUnit);
  fSetP1Cmd->SetParameter(p1xValFcn);
  fSetP1Cmd->SetParameter(p1xValBinScheme);
  fSetP1Cmd->SetParameter(p1yValMin);
  fSetP1Cmd->SetParameter(p1yValMax);
  fSetP1Cmd->SetParameter(p1yValUnit);
  fSetP1Cmd->SetParameter(p1yValFcn);
  fSetP1Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  

//
// public functions
//

//_____________________________________________________________________________
void G4P1Messenger::SetNewValue(G4UIcommand* command, G4String newValues)
{
  // tokenize parameters in a vector
  std::vector<G4String> parameters;
  G4Analysis::Tokenize(newValues, parameters);
  // check consistency
  if ( G4int(parameters.size()) != command->GetParameterEntries() ) {
    // Should never happen but let's check anyway for consistency
    fHelper->WarnAboutParameters(command, parameters.size());
    return;
  }  

  if ( command == fCreateP1Cmd.get() ) {  
    auto counter = 0;
    auto name = parameters[counter++];
    auto title = parameters[counter++];
    G4AnalysisMessengerHelper::BinData xdata;
    fHelper->GetBinData(xdata, parameters, counter);
    auto xunit = GetUnitValue(xdata.fSunit);
    G4AnalysisMessengerHelper::ValueData ydata;
    fHelper->GetValueData(ydata, parameters, counter);
    auto yunit = GetUnitValue(ydata.fSunit);
    fManager->CreateP1(name, title, 
                       xdata.fNbins, xdata.fVmin*xunit, xdata.fVmax*xunit, 
                       ydata.fVmin*yunit, ydata.fVmax*yunit, 
                       xdata.fSunit, ydata.fSunit, 
                       xdata.fSfcn, ydata.fSfcn, 
                       xdata.fSbinScheme);     
  }
  else if ( command == fSetP1Cmd.get() ) { 
    auto counter = 0;
    auto id = G4UIcommand::ConvertToInt(parameters[counter++]);
    G4AnalysisMessengerHelper::BinData xdata;
    fHelper->GetBinData(xdata, parameters, counter);
    auto xunit = GetUnitValue(xdata.fSunit);
    G4AnalysisMessengerHelper::ValueData ydata;
    fHelper->GetValueData(ydata, parameters, counter);
    auto yunit = GetUnitValue(ydata.fSunit);
    fManager->SetP1(id, 
                    xdata.fNbins, xdata.fVmin*xunit, xdata.fVmax*xunit, 
                    ydata.fVmin*yunit, ydata.fVmax*yunit, 
                    xdata.fSunit, ydata.fSunit, 
                    xdata.fSfcn, ydata.fSfcn, 
                    xdata.fSbinScheme);     
  }
  else if ( command == fSetP1XCmd.get() ) { 
    // Save values
    auto counter = 0;
    fXId = G4UIcommand::ConvertToInt(parameters[counter++]);
    fHelper->GetBinData(fXData, parameters, counter);
    // Set values
    // (another set may follow if setY is also called)
    auto xunit = GetUnitValue(fXData.fSunit);
    fManager->SetP1(fXId, 
                    fXData.fNbins, fXData.fVmin*xunit, fXData.fVmax*xunit, 
                    0., 0., 
                    fXData.fSunit, "none", 
                    fXData.fSfcn, "none", 
                    fXData.fSbinScheme);     
  }
  else if ( command == fSetP1YCmd.get() ) { 
    // Check if setX command was called
    auto counter = 0;
    auto id = G4UIcommand::ConvertToInt(parameters[counter++]);
    if ( fXId == -1 || fXId != id ) {
      fHelper->WarnAboutSetCommands();
      return;
    }
    auto xunit = GetUnitValue(fXData.fSunit);
    G4AnalysisMessengerHelper::ValueData ydata;
    fHelper->GetValueData(ydata, parameters, counter);
    auto yunit = GetUnitValue(ydata.fSunit);
    fManager->SetP1(id, 
                    fXData.fNbins, fXData.fVmin*xunit, fXData.fVmax*xunit,
                    ydata.fVmin*yunit, ydata.fVmax*yunit, 
                    fXData.fSunit, ydata.fSunit, 
                    fXData.fSfcn, ydata.fSfcn, 
                    fXData.fSbinScheme);     
    fXId = -1;                        
  }
  else if ( command == fSetP1TitleCmd.get() ) { 
    auto counter = 0;
    auto id = G4UIcommand::ConvertToInt(parameters[counter++]); 
    auto title = parameters[counter++];
    fManager->SetP1Title(id, title);     
  }
  else if ( command == fSetP1XAxisCmd.get() ) { 
    auto counter = 0;
    auto id = G4UIcommand::ConvertToInt(parameters[counter++]); 
    auto xaxis = parameters[counter++];
    fManager->SetP1XAxisTitle(id, xaxis);     
  }
  else if ( command == fSetP1YAxisCmd.get() ) { 
    auto counter = 0;
    auto id = G4UIcommand::ConvertToInt(parameters[counter++]); 
    auto yaxis = parameters[counter++];
    fManager->SetP1YAxisTitle(id, yaxis);     
  }
}  
