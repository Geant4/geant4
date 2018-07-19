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

#include "G4P2Messenger.hh"
#include "G4VAnalysisManager.hh"
#include "G4AnalysisUtilities.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"

#include <iostream>

using namespace G4Analysis;

//_____________________________________________________________________________
G4P2Messenger::G4P2Messenger(G4VAnalysisManager* manager)
  : G4UImessenger(),
    fManager(manager),
    fHelper(nullptr),
    fDirectory(nullptr),
    fCreateP2Cmd(nullptr),
    fSetP2Cmd(nullptr),
    fSetP2XCmd(nullptr),
    fSetP2YCmd(nullptr),
    fSetP2ZCmd(nullptr),
    fSetP2TitleCmd(nullptr), 
    fSetP2XAxisCmd(nullptr), 
    fSetP2YAxisCmd(nullptr)
{  
  fHelper = G4Analysis::make_unique<G4AnalysisMessengerHelper>("p2");

  fDirectory = fHelper->CreateHnDirectory();

  CreateP2Cmd();

  SetP2Cmd();
  fSetP2XCmd = fHelper->CreateSetBinsCommand("x", this);
  fSetP2YCmd = fHelper->CreateSetBinsCommand("y", this);
  fSetP2ZCmd = fHelper->CreateSetValuesCommand("z", this);

  fSetP2TitleCmd = fHelper->CreateSetTitleCommand(this);
  fSetP2XAxisCmd = fHelper->CreateSetAxisCommand("x", this);
  fSetP2YAxisCmd = fHelper->CreateSetAxisCommand("y", this);
  fSetP2ZAxisCmd = fHelper->CreateSetAxisCommand("z", this);
}

//_____________________________________________________________________________
G4P2Messenger::~G4P2Messenger()
{}

//
// private functions
//

//_____________________________________________________________________________
void G4P2Messenger::CreateP2Cmd()
{
  auto p2Name = new G4UIparameter("name", 's', false);
  p2Name->SetGuidance("Profile name (label)");
  
  auto p2Title = new G4UIparameter("title", 's', false);
  p2Title->SetGuidance("Profile title");

  auto p2xNbins0 = new G4UIparameter("xnbins0", 'i', true);
  p2xNbins0->SetGuidance("Number of x-bins (default = 100)");
  p2xNbins0->SetGuidance("Can be reset with /analysis/p2/set command");
  p2xNbins0->SetDefaultValue(100);
  
  auto p2xValMin0 = new G4UIparameter("xvalMin0", 'd', true);
  p2xValMin0->SetGuidance("Minimum x-value, expressed in unit (default = 0.)");
  p2xValMin0->SetGuidance("Can be reset with /analysis/p2/set command");
  p2xValMin0->SetDefaultValue(0.);
  
  auto p2xValMax0 = new G4UIparameter("xvalMax0", 'd', true);
  p2xValMax0->SetGuidance("Maximum x-value, expressed in unit (default = 1.)");
  p2xValMax0->SetGuidance("Can be reset with /analysis/p2/set command");
  p2xValMax0->SetDefaultValue(1.);

  auto p2xValUnit0 = new G4UIparameter("xvalUnit0", 's', true);
  p2xValUnit0->SetGuidance("The unit applied to filled x-values and xvalMin0, xvalMax0");
  p2xValUnit0->SetDefaultValue("none");
  
  auto p2xValFcn0 = new G4UIparameter("xvalFcn0", 's', true);
  G4String fcnxGuidance = "The function applied to filled x-values (log, log10, exp, none).";
  p2xValFcn0->SetGuidance(fcnxGuidance);
  p2xValFcn0->SetParameterCandidates("log log10 exp none");
  p2xValFcn0->SetDefaultValue("none");
    
  auto p2xValBinScheme0 = new G4UIparameter("xvalBinScheme0", 's', true);
  G4String binSchemeGuidance = "The binning scheme (linear, log).";
  p2xValBinScheme0->SetParameterCandidates("linear log");
  p2xValBinScheme0->SetGuidance(binSchemeGuidance);
  p2xValBinScheme0->SetDefaultValue("linear");
  
  auto p2yNbins0 = new G4UIparameter("ynbins0", 'i', true);
  p2yNbins0->SetGuidance("Number of y-bins (default = 100)");
  p2yNbins0->SetGuidance("Can be reset with /analysis/p2/set command");
  p2yNbins0->SetDefaultValue(100);
  
  auto p2yValMin0 = new G4UIparameter("yvalMin0", 'd', true);
  p2yValMin0->SetGuidance("Minimum y-value, expressed in unit (default = 0.)");
  p2yValMin0->SetGuidance("Can be reset with /analysis/p2/set command");
  p2yValMin0->SetDefaultValue(0.);
  
  auto p2yValMax0 = new G4UIparameter("yvalMax0", 'd', true);
  p2yValMax0->SetGuidance("Maximum y-value, expressed in unit (default = 1.)");
  p2yValMax0->SetGuidance("Can be reset with /analysis/p2/set command");
  p2yValMax0->SetDefaultValue(1.);

  auto p2yValUnit0 = new G4UIparameter("yvalUnit0", 's', true);
  p2yValUnit0->SetGuidance("The unit applied to filled y-values and yvalMin0, yvalMax0");
  p2yValUnit0->SetDefaultValue("none");
  
  auto p2yValFcn0 = new G4UIparameter("yvalFcn0", 's', true);
  G4String fcnyGuidance = "The function applied to filled y-values (log, log10, exp, none).";
  p2yValFcn0->SetGuidance(fcnyGuidance);
  p2yValFcn0->SetParameterCandidates("log log10 exp none");
  p2yValFcn0->SetDefaultValue("none");
    
  auto p2yValBinScheme0 = new G4UIparameter("yvalBinScheme0", 's', true);
  p2yValBinScheme0->SetParameterCandidates("linear log");
  p2yValBinScheme0->SetGuidance(binSchemeGuidance);
  p2yValBinScheme0->SetDefaultValue("linear");
  
  auto p2zValMin0 = new G4UIparameter("zvalMin0", 'd', true);
  p2zValMin0->SetGuidance("Minimum z-value, expressed in unit (default = 0.)");
  p2zValMin0->SetGuidance("Can be reset with /analysis/p2/set command");
  p2zValMin0->SetDefaultValue(0.);
  
  auto p2zValMax0 = new G4UIparameter("zvalMax0", 'd', true);
  p2zValMax0->SetGuidance("Maximum z-value, expressed in unit (default = 1.)");
  p2zValMax0->SetGuidance("Can be reset with /analysis/p2/set command");
  p2zValMax0->SetDefaultValue(1.);

  auto p2zValUnit0 = new G4UIparameter("zvalUnit0", 's', true);
  p2zValUnit0->SetGuidance("The unit applied to filled z-values and zvalMin0, zvalMax0");
  p2zValUnit0->SetDefaultValue("none");
  
  auto p2zValFcn0 = new G4UIparameter("zvalFcn0", 's', true);
  G4String fcnzGuidance = "The function applied to filled z-values (log, log10, exp, none).";
  p2zValFcn0->SetGuidance(fcnzGuidance);
  p2zValFcn0->SetParameterCandidates("log log10 exp none");
  p2zValFcn0->SetDefaultValue("none");
  
  fCreateP2Cmd = G4Analysis::make_unique<G4UIcommand>("/analysis/p2/create", this);
  fCreateP2Cmd->SetGuidance("Create 2D profile");
  fCreateP2Cmd->SetParameter(p2Name);
  fCreateP2Cmd->SetParameter(p2Title);
  fCreateP2Cmd->SetParameter(p2xNbins0);
  fCreateP2Cmd->SetParameter(p2xValMin0);
  fCreateP2Cmd->SetParameter(p2xValMax0);
  fCreateP2Cmd->SetParameter(p2xValUnit0);
  fCreateP2Cmd->SetParameter(p2xValFcn0);
  fCreateP2Cmd->SetParameter(p2xValBinScheme0);
  fCreateP2Cmd->SetParameter(p2yNbins0);
  fCreateP2Cmd->SetParameter(p2yValMin0);
  fCreateP2Cmd->SetParameter(p2yValMax0);
  fCreateP2Cmd->SetParameter(p2yValUnit0);
  fCreateP2Cmd->SetParameter(p2yValFcn0);
  fCreateP2Cmd->SetParameter(p2yValBinScheme0);
  fCreateP2Cmd->SetParameter(p2zValMin0);
  fCreateP2Cmd->SetParameter(p2zValMax0);
  fCreateP2Cmd->SetParameter(p2zValUnit0);
  fCreateP2Cmd->SetParameter(p2zValFcn0);
  fCreateP2Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  


//_____________________________________________________________________________
void G4P2Messenger::SetP2Cmd()
{
  auto p2Id = new G4UIparameter("id", 'i', false);
  p2Id->SetGuidance("Profile id");
  p2Id->SetParameterRange("id>=0");
  
  auto p2xNbins = new G4UIparameter("xnbins", 'i', false);
  p2xNbins->SetGuidance("Number of x-bins");
  
  auto p2xValMin = new G4UIparameter("xvalMin", 'd', false);
  p2xValMin->SetGuidance("Minimum x-value, expressed in unit");
  
  auto p2xValMax = new G4UIparameter("xvalMax", 'd', false);
  p2xValMax->SetGuidance("Maximum x-value, expressed in unit");
  
  auto p2xValUnit = new G4UIparameter("xvalUnit", 's', true);
  p2xValUnit->SetGuidance("The unit applied to filled x-values and xvalMin0, xvalMax0");
  p2xValUnit->SetDefaultValue("none");
 
  auto p2xValFcn = new G4UIparameter("xvalFcn", 's', true);
  p2xValFcn->SetParameterCandidates("log log10 exp none");
  G4String fcnxGuidance = "The function applied to filled x-values (log, log10, exp, none).";
  p2xValFcn->SetGuidance(fcnxGuidance);
  p2xValFcn->SetDefaultValue("none");
    
  auto p2xValBinScheme = new G4UIparameter("xvalBinScheme", 's', true);
  G4String binSchemeGuidance = "The binning scheme (linear, log).";
  p2xValBinScheme->SetParameterCandidates("linear log");
  p2xValBinScheme->SetGuidance(binSchemeGuidance);
  p2xValBinScheme->SetDefaultValue("linear");
 
  auto p2yNbins = new G4UIparameter("nybins", 'i', true);
  p2yNbins->SetGuidance("Number of y-bins");
  
  auto p2yValMin = new G4UIparameter("yvalMin", 'd', true);
  p2yValMin->SetGuidance("Minimum y-value, expressed in unit");
  
  auto p2yValMax = new G4UIparameter("yvalMax", 'd', true);
  p2yValMax->SetGuidance("Maximum y-value, expressed in unit");
  
  auto p2yValUnit = new G4UIparameter("yvalUnit", 's', true);
  p2yValUnit->SetGuidance("The unit applied to filled y-values and yvalMin0, yvalMax0");
  p2yValUnit->SetDefaultValue("none");
 
  auto p2yValFcn = new G4UIparameter("yvalFcn", 's', true);
  p2yValFcn->SetParameterCandidates("log log10 exp none");
  G4String fcnyGuidance = "The function applied to filled y-values (log, log10, exp, none).";
  p2yValFcn->SetGuidance(fcnyGuidance);
  p2yValFcn->SetDefaultValue("none");
    
  auto p2yValBinScheme = new G4UIparameter("yvalBinScheme", 's', true);
  p2yValBinScheme->SetParameterCandidates("linear log");
  p2yValBinScheme->SetGuidance(binSchemeGuidance);
  p2yValBinScheme->SetDefaultValue("linear");
 
  auto p2zValMin = new G4UIparameter("zvalMin", 'd', true);
  p2zValMin->SetGuidance("Minimum z-value, expressed in unit");
  
  auto p2zValMax = new G4UIparameter("zvalMax", 'd', true);
  p2zValMax->SetGuidance("Maximum z-value, expressed in unit");
  
  auto p2zValUnit = new G4UIparameter("zvalUnit", 's', true);
  p2zValUnit->SetGuidance("The unit applied to filled z-values and zvalMin0, zvalMax0");
  p2zValUnit->SetDefaultValue("none");
 
  auto p2zValFcn = new G4UIparameter("zvalFcn", 's', true);
  p2zValFcn->SetParameterCandidates("log log10 exp none");
  G4String fcnzGuidance = "The function applied to filled z-values (log, log10, exp, none).";
  p2zValFcn->SetGuidance(fcnzGuidance);
  p2zValFcn->SetDefaultValue("none");
 
  fSetP2Cmd = G4Analysis::make_unique<G4UIcommand>("/analysis/p2/set", this);
  fSetP2Cmd->SetGuidance("Set parameters for the 2D profile of given id:");
  fSetP2Cmd->SetGuidance("  nxbins; xvalMin; xvalMax; xunit; xbinScheme");
  fSetP2Cmd->SetGuidance("  nybins; yvalMin; yvalMax; yunit; ybinScheme");
  fSetP2Cmd->SetGuidance("  zvalMin; zvalMax; zunit; zfunction");
  fSetP2Cmd->SetParameter(p2Id);
  fSetP2Cmd->SetParameter(p2xNbins);
  fSetP2Cmd->SetParameter(p2xValMin);
  fSetP2Cmd->SetParameter(p2xValMax);
  fSetP2Cmd->SetParameter(p2xValUnit);
  fSetP2Cmd->SetParameter(p2xValFcn);
  fSetP2Cmd->SetParameter(p2xValBinScheme);
  fSetP2Cmd->SetParameter(p2yNbins);
  fSetP2Cmd->SetParameter(p2yValMin);
  fSetP2Cmd->SetParameter(p2yValMax);
  fSetP2Cmd->SetParameter(p2yValUnit);
  fSetP2Cmd->SetParameter(p2yValFcn);
  fSetP2Cmd->SetParameter(p2yValBinScheme);
  fSetP2Cmd->SetParameter(p2zValMin);
  fSetP2Cmd->SetParameter(p2zValMax);
  fSetP2Cmd->SetParameter(p2zValUnit);
  fSetP2Cmd->SetParameter(p2zValFcn);
  fSetP2Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  


//
// public functions
//

//_____________________________________________________________________________
void G4P2Messenger::SetNewValue(G4UIcommand* command, G4String newValues)
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

  if ( command == fCreateP2Cmd.get() ) {  
    auto counter = 0;
    auto name = parameters[counter++];
    auto title = parameters[counter++];
    G4AnalysisMessengerHelper::BinData xdata;
    fHelper->GetBinData(xdata, parameters, counter);
    auto xunit = GetUnitValue(xdata.fSunit);
    G4AnalysisMessengerHelper::BinData ydata;
    fHelper->GetBinData(ydata, parameters, counter);
    auto yunit = GetUnitValue(ydata.fSunit);
    G4AnalysisMessengerHelper::ValueData zdata;
    fHelper->GetValueData(zdata, parameters, counter);
    auto zunit = GetUnitValue(zdata.fSunit);
    fManager->CreateP2(name, title, 
                       xdata.fNbins, xdata.fVmin*xunit, xdata.fVmax*xunit, 
                       ydata.fNbins, ydata.fVmin*yunit, ydata.fVmax*yunit, 
                       zdata.fVmin*zunit, zdata.fVmax*zunit, 
                       xdata.fSunit, ydata.fSunit, zdata.fSunit,
                       xdata.fSfcn, ydata.fSfcn, zdata.fSfcn,
                       xdata.fSbinScheme, ydata.fSbinScheme);     
  }
  else if ( command == fSetP2Cmd.get() ) { 
    auto counter = 0;
    auto id = G4UIcommand::ConvertToInt(parameters[counter++]);
    G4AnalysisMessengerHelper::BinData xdata;
    fHelper->GetBinData(xdata, parameters, counter);
    auto xunit = GetUnitValue(xdata.fSunit);
    G4AnalysisMessengerHelper::BinData ydata;
    fHelper->GetBinData(ydata, parameters, counter);
    auto yunit = GetUnitValue(ydata.fSunit);
    G4AnalysisMessengerHelper::ValueData zdata;
    fHelper->GetValueData(zdata, parameters, counter);
    auto zunit = GetUnitValue(zdata.fSunit);
    fManager->SetP2(id, 
                    xdata.fNbins, xdata.fVmin*xunit, xdata.fVmax*xunit, 
                    ydata.fNbins, ydata.fVmin*yunit, ydata.fVmax*yunit, 
                    zdata.fVmin*zunit, zdata.fVmax*zunit, 
                    xdata.fSunit, ydata.fSunit, zdata.fSunit,
                    xdata.fSfcn, ydata.fSfcn, zdata.fSfcn,
                    xdata.fSbinScheme, ydata.fSbinScheme);     
  }
  else if ( command == fSetP2XCmd.get() ) { 
    // Only save values
    auto counter = 0;
    fXId = G4UIcommand::ConvertToInt(parameters[counter++]);
    fHelper->GetBinData(fXData, parameters, counter);
  }
  else if ( command == fSetP2YCmd.get() ) { 
    // Save values
    auto counter = 0;
    fYId = G4UIcommand::ConvertToInt(parameters[counter++]);
    // Check if setX command was called
    if ( fXId == -1 || fXId != fYId ) {
      fHelper->WarnAboutSetCommands();
      return;
    }
    fHelper->GetBinData(fYData, parameters, counter);
    // Set values
    // (another set may follow if setZ is also called)
    auto xunit = GetUnitValue(fXData.fSunit);
    auto yunit = GetUnitValue(fYData.fSunit);
    fManager->SetP2(fYId, 
                    fXData.fNbins, fXData.fVmin*xunit, fXData.fVmax*xunit,
                    fYData.fNbins, fYData.fVmin*yunit, fYData.fVmax*yunit,
                    0., 0., 
                    fXData.fSunit, fYData.fSunit, "none",
                    fXData.fSfcn, fYData.fSfcn, "none",
                    fXData.fSbinScheme, fYData.fSbinScheme);     
  }
  else if ( command == fSetP2ZCmd.get() ) { 
    auto counter = 0;
    auto id = G4UIcommand::ConvertToInt(parameters[counter++]);
    // Check if setX and setY command was called
    if ( fXId == -1 || fXId != id || 
         fYId == -1 || fYId != id ) {
      fHelper->WarnAboutSetCommands();
      return;
    }
    auto xunit = GetUnitValue(fXData.fSunit);
    auto yunit = GetUnitValue(fYData.fSunit);
    G4AnalysisMessengerHelper::ValueData zdata;
    fHelper->GetValueData(zdata, parameters, counter);
    auto zunit = GetUnitValue(zdata.fSunit);
    fManager->SetP2(id, 
                    fXData.fNbins, fXData.fVmin*xunit, fXData.fVmax*xunit,
                    fYData.fNbins, fYData.fVmin*yunit, fYData.fVmax*yunit,
                    zdata.fVmin*zunit, zdata.fVmax*zunit, 
                    fXData.fSunit, fYData.fSunit, zdata.fSunit,
                    fXData.fSfcn, fYData.fSfcn, zdata.fSfcn,
                    fXData.fSbinScheme, fYData.fSbinScheme);     
    fXId = -1;                        
    fYId = -1;                        
  }
  else if ( command == fSetP2TitleCmd.get() ) { 
    auto counter = 0;
    auto id = G4UIcommand::ConvertToInt(parameters[counter++]); 
    auto title = parameters[counter++];
    fManager->SetP2Title(id, title);     
  }
  else if ( command == fSetP2XAxisCmd.get() ) { 
    auto counter = 0;
    auto id = G4UIcommand::ConvertToInt(parameters[counter++]); 
    auto xaxis = parameters[counter++];
    fManager->SetP2XAxisTitle(id, xaxis);     
  }
  else if ( command == fSetP2YAxisCmd.get() ) { 
    auto counter = 0;
    auto id = G4UIcommand::ConvertToInt(parameters[counter++]); 
    auto yaxis = parameters[counter++];
    fManager->SetP2YAxisTitle(id, yaxis);     
  }
  else if ( command == fSetP2ZAxisCmd.get() ) { 
    auto counter = 0;
    auto id = G4UIcommand::ConvertToInt(parameters[counter++]); 
    auto zaxis = parameters[counter++];
    fManager->SetP2ZAxisTitle(id, zaxis);     
  }
}  
