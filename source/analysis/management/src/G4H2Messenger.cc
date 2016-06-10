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
// $Id: G4H2Messenger.cc 66310 2012-12-17 11:56:35Z ihrivnac $

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#include "G4H2Messenger.hh"
#include "G4VAnalysisManager.hh"
#include "G4AnalysisUtilities.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"

#include <iostream>

using namespace G4Analysis;

//_____________________________________________________________________________
G4H2Messenger::G4H2Messenger(G4VAnalysisManager* manager)
  : G4UImessenger(),
    fManager(manager),
    fHelper(nullptr),
    fDirectory(nullptr),
    fCreateH2Cmd(nullptr),
    fSetH2Cmd(nullptr),
    fSetH2XCmd(nullptr),
    fSetH2YCmd(nullptr),
    fSetH2TitleCmd(nullptr), 
    fSetH2XAxisCmd(nullptr), 
    fSetH2YAxisCmd(nullptr),
    fXId(-1),
    fXData()
{  
  fHelper = G4Analysis::make_unique<G4AnalysisMessengerHelper>("h2");

  fDirectory = fHelper->CreateHnDirectory();

  CreateH2Cmd();

  SetH2Cmd();
  fSetH2XCmd = fHelper->CreateSetBinsCommand("x", this);
  fSetH2YCmd = fHelper->CreateSetBinsCommand("y", this);
  
  fSetH2TitleCmd = fHelper->CreateSetTitleCommand(this);
  fSetH2XAxisCmd = fHelper->CreateSetAxisCommand("x", this);
  fSetH2YAxisCmd = fHelper->CreateSetAxisCommand("y", this);
  fSetH2ZAxisCmd = fHelper->CreateSetAxisCommand("z", this);
}

//_____________________________________________________________________________
G4H2Messenger::~G4H2Messenger()
{}

//
// private functions
//

//_____________________________________________________________________________
void G4H2Messenger::CreateH2Cmd()
{
  auto h2Name = new G4UIparameter("name", 's', false);
  h2Name->SetGuidance("Histogram name (label)");
  
  auto h2Title = new G4UIparameter("title", 's', false);
  h2Title->SetGuidance("Histogram title");

  auto h2xNbins0 = new G4UIparameter("xnbins0", 'i', true);
  h2xNbins0->SetGuidance("Number of x-bins (default = 100)");
  h2xNbins0->SetGuidance("Can be reset with /analysis/h2/set command");
  h2xNbins0->SetDefaultValue(100);
  
  auto h2xValMin0 = new G4UIparameter("xvalMin0", 'd', true);
  h2xValMin0->SetGuidance("Minimum x-value, expressed in unit (default = 0.)");
  h2xValMin0->SetGuidance("Can be reset with /analysis/h2/set command");
  h2xValMin0->SetDefaultValue(0.);
  
  auto h2xValMax0 = new G4UIparameter("xvalMax0", 'd', true);
  h2xValMax0->SetGuidance("Maximum x-value, expressed in unit (default = 1.)");
  h2xValMax0->SetGuidance("Can be reset with /analysis/h2/set command");
  h2xValMax0->SetDefaultValue(1.);

  auto h2xValUnit0 = new G4UIparameter("xvalUnit0", 's', true);
  h2xValUnit0->SetGuidance("The unit applied to filled x-values and xvalMin0, xvalMax0");
  h2xValUnit0->SetDefaultValue("none");
  
  auto h2xValFcn0 = new G4UIparameter("xvalFcn0", 's', true);
  G4String fcnxGuidance = "The function applied to filled x-values (log, log10, exp, none).";
  h2xValFcn0->SetGuidance(fcnxGuidance);
  h2xValFcn0->SetParameterCandidates("log log10 exp none");
  h2xValFcn0->SetDefaultValue("none");
  
  auto h2xValBinScheme0 = new G4UIparameter("xvalBinScheme0", 's', true);
  G4String xbinSchemeGuidance = "The binning scheme (linear, log).";
  h2xValBinScheme0->SetParameterCandidates("linear log");
  h2xValBinScheme0->SetGuidance(xbinSchemeGuidance);
  h2xValBinScheme0->SetDefaultValue("linear");
  
  auto h2yNbins0 = new G4UIparameter("ynbins0", 'i', true);
  h2yNbins0->SetGuidance("Number of y-bins (default = 100)");
  h2yNbins0->SetGuidance("Can be reset with /analysis/h2/set command");
  h2yNbins0->SetDefaultValue(100);
  
  auto h2yValMin0 = new G4UIparameter("yvalMin0", 'd', true);
  h2yValMin0->SetGuidance("Minimum y-value, expressed in unit (default = 0.)");
  h2yValMin0->SetGuidance("Can be reset with /analysis/h2/set command");
  h2yValMin0->SetDefaultValue(0.);
  
  auto h2yValMax0 = new G4UIparameter("yvalMax0", 'd', true);
  h2yValMax0->SetGuidance("Maximum y-value, expressed in unit (default = 1.)");
  h2yValMax0->SetGuidance("Can be reset with /analysis/h2/set command");
  h2yValMax0->SetDefaultValue(1.);

  auto h2yValUnit0 = new G4UIparameter("yvalUnit0", 's', true);
  h2yValUnit0->SetGuidance("The unit applied to filled y-values and yvalMin0, yvalMax0");
  h2yValUnit0->SetDefaultValue("none");
  
  auto h2yValFcn0 = new G4UIparameter("yvalFcn0", 's', true);
  G4String fcnyGuidance = "The function applied to filled y-values (log, log10, exp, none).";
  h2yValFcn0->SetGuidance(fcnyGuidance);
  h2yValFcn0->SetParameterCandidates("log log10 exp none");
  h2yValFcn0->SetDefaultValue("none");

  auto h2yValBinScheme0 = new G4UIparameter("yvalBinScheme0", 's', true);
  G4String ybinSchemeGuidance = "The binning scheme (linear, log).";
  h2yValBinScheme0->SetParameterCandidates("linear log");
  h2yValBinScheme0->SetGuidance(ybinSchemeGuidance);
  h2yValBinScheme0->SetDefaultValue("linear");
  
  fCreateH2Cmd = G4Analysis::make_unique<G4UIcommand>("/analysis/h2/create", this);
  fCreateH2Cmd->SetGuidance("Create 2D histogram");
  fCreateH2Cmd->SetParameter(h2Name);
  fCreateH2Cmd->SetParameter(h2Title);
  fCreateH2Cmd->SetParameter(h2xNbins0);
  fCreateH2Cmd->SetParameter(h2xValMin0);
  fCreateH2Cmd->SetParameter(h2xValMax0);
  fCreateH2Cmd->SetParameter(h2xValUnit0);
  fCreateH2Cmd->SetParameter(h2xValFcn0);
  fCreateH2Cmd->SetParameter(h2xValBinScheme0);
  fCreateH2Cmd->SetParameter(h2yNbins0);
  fCreateH2Cmd->SetParameter(h2yValMin0);
  fCreateH2Cmd->SetParameter(h2yValMax0);
  fCreateH2Cmd->SetParameter(h2yValUnit0);
  fCreateH2Cmd->SetParameter(h2yValFcn0);
  fCreateH2Cmd->SetParameter(h2yValBinScheme0);
  fCreateH2Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  


//_____________________________________________________________________________
void G4H2Messenger::SetH2Cmd()
{
  auto h2Id = new G4UIparameter("id", 'i', false);
  h2Id->SetGuidance("Histogram id");
  h2Id->SetParameterRange("id>=0");
  
  auto h2xNbins = new G4UIparameter("xnbins", 'i', false);
  h2xNbins->SetGuidance("Number of x-bins");
  
  auto h2xValMin = new G4UIparameter("xvalMin", 'd', false);
  h2xValMin->SetGuidance("Minimum x-value, expressed in unit");
  
  auto h2xValMax = new G4UIparameter("xvalMax", 'd', false);
  h2xValMax->SetGuidance("Maximum x-value, expressed in unit");
  
  auto h2xValUnit = new G4UIparameter("xvalUnit", 's', false);
  h2xValUnit->SetGuidance("The unit applied to filled x-values and xvalMin, xvalMax");
  h2xValUnit->SetDefaultValue("none");
 
  auto h2xValFcn = new G4UIparameter("xvalFcn", 's', false);
  h2xValFcn->SetParameterCandidates("log log10 exp none");
  G4String fcnxGuidance = "The function applied to filled x-values (log, log10, exp, none).";
  h2xValFcn->SetGuidance(fcnxGuidance);
  h2xValFcn->SetDefaultValue("none");
 
  auto h2xValBinScheme = new G4UIparameter("xvalBinScheme", 's', true);
  G4String xbinSchemeGuidance = "The binning scheme (linear, log).";
  h2xValBinScheme->SetParameterCandidates("linear log");
  h2xValBinScheme->SetGuidance(xbinSchemeGuidance);
  h2xValBinScheme->SetDefaultValue("linear");
  
  auto h2yNbins = new G4UIparameter("nybins", 'i', false);
  h2yNbins->SetGuidance("Number of y-bins");
  
  auto h2yValMin = new G4UIparameter("yvalMin", 'd', false);
  h2yValMin->SetGuidance("Minimum y-value, expressed in unit");
  
  auto h2yValMax = new G4UIparameter("yvalMax", 'd', false);
  h2yValMax->SetGuidance("Maximum y-value, expressed in unit");
  
  auto h2yValUnit = new G4UIparameter("yvalUnit", 's', true);
  h2yValUnit->SetGuidance("The unit applied to filled y-values and yvalMin, yvalMax");
  h2yValUnit->SetDefaultValue("none");
 
  auto h2yValFcn = new G4UIparameter("yvalFcn", 's', false);
  h2yValFcn->SetParameterCandidates("log log10 exp none");
  G4String fcnyGuidance = "The function applied to filled y-values (log, log10, exp, none).";
  h2yValFcn->SetGuidance(fcnyGuidance);
  h2yValFcn->SetDefaultValue("none");
 
  auto h2yValBinScheme = new G4UIparameter("yvalBinScheme", 's', true);
  G4String ybinSchemeGuidance = "The binning scheme (linear, log).";
  h2yValBinScheme->SetParameterCandidates("linear log");
  h2yValBinScheme->SetGuidance(ybinSchemeGuidance);
  h2yValBinScheme->SetDefaultValue("linear");

  fSetH2Cmd = G4Analysis::make_unique<G4UIcommand>("/analysis/h2/set", this);
  fSetH2Cmd->SetGuidance("Set parameters for the 2D histogram of given id:");
  fSetH2Cmd->SetGuidance("  nxbins; xvalMin; xvalMax; xunit; xfunction; xbinScheme");
  fSetH2Cmd->SetGuidance("  nybins; yvalMin; yvalMax; yunit; yfunction; ybinScheme");
  fSetH2Cmd->SetParameter(h2Id);
  fSetH2Cmd->SetParameter(h2xNbins);
  fSetH2Cmd->SetParameter(h2xValMin);
  fSetH2Cmd->SetParameter(h2xValMax);
  fSetH2Cmd->SetParameter(h2xValUnit);
  fSetH2Cmd->SetParameter(h2xValFcn);
  fSetH2Cmd->SetParameter(h2xValBinScheme);
  fSetH2Cmd->SetParameter(h2yNbins);
  fSetH2Cmd->SetParameter(h2yValMin);
  fSetH2Cmd->SetParameter(h2yValMax);
  fSetH2Cmd->SetParameter(h2yValUnit);
  fSetH2Cmd->SetParameter(h2yValFcn);
  fSetH2Cmd->SetParameter(h2yValBinScheme);
  fSetH2Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  


//
// public functions
//

//_____________________________________________________________________________
void G4H2Messenger::SetNewValue(G4UIcommand* command, G4String newValues)
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

  if ( command == fCreateH2Cmd.get() ) {  
    auto counter = 0;
    auto name = parameters[counter++];
    auto title = parameters[counter++];
    G4AnalysisMessengerHelper::BinData xdata;
    fHelper->GetBinData(xdata, parameters, counter);
    auto xunit = GetUnitValue(xdata.fSunit);
    G4AnalysisMessengerHelper::BinData ydata;
    fHelper->GetBinData(ydata, parameters, counter);
    auto yunit = GetUnitValue(ydata.fSunit);
    fManager->CreateH2(name, title, 
                       xdata.fNbins, xdata.fVmin*xunit, xdata.fVmax*xunit, 
                       ydata.fNbins, ydata.fVmin*yunit, ydata.fVmax*yunit, 
                       xdata.fSunit, ydata.fSunit, 
                       xdata.fSfcn, ydata.fSfcn, 
                       xdata.fSbinScheme, ydata.fSbinScheme);     
  }
  else if ( command == fSetH2Cmd.get() ) { 
    auto counter = 0;
    auto id = G4UIcommand::ConvertToInt(parameters[counter++]);
    G4AnalysisMessengerHelper::BinData xdata;
    fHelper->GetBinData(xdata, parameters, counter);
    auto xunit = GetUnitValue(xdata.fSunit);
    G4AnalysisMessengerHelper::BinData ydata;
    fHelper->GetBinData(ydata, parameters, counter);
    auto yunit = GetUnitValue(ydata.fSunit);
    fManager->SetH2(id, 
                    xdata.fNbins, xdata.fVmin*xunit, xdata.fVmax*xunit, 
                    ydata.fNbins, ydata.fVmin*yunit, ydata.fVmax*yunit, 
                    xdata.fSunit, ydata.fSunit, 
                    xdata.fSfcn, ydata.fSfcn, 
                    xdata.fSbinScheme, ydata.fSbinScheme);     
  }
  else if ( command == fSetH2XCmd.get() ) { 
    // Only save values
    auto counter = 0;
    fXId = G4UIcommand::ConvertToInt(parameters[counter++]);
    fHelper->GetBinData(fXData, parameters, counter);
  }
  else if ( command == fSetH2YCmd.get() ) { 
    // Check if setX command was called
    auto counter = 0;
    auto id = G4UIcommand::ConvertToInt(parameters[counter++]);
    if ( fXId == -1 || fXId != id ) {
      fHelper->WarnAboutSetCommands();
      return;
    }
    auto xunit = GetUnitValue(fXData.fSunit);
    G4AnalysisMessengerHelper::BinData ydata;
    fHelper->GetBinData(ydata, parameters, counter);
    auto yunit = GetUnitValue(ydata.fSunit);
    fManager->SetH2(id, 
                    fXData.fNbins, fXData.fVmin*xunit, fXData.fVmax*xunit,
                    ydata.fNbins, ydata.fVmin*yunit, ydata.fVmax*yunit, 
                    fXData.fSunit, ydata.fSunit, 
                    fXData.fSfcn, ydata.fSfcn,  
                    fXData.fSbinScheme, ydata.fSbinScheme); 
    fXId = -1;                        
  }
  else if ( command == fSetH2TitleCmd.get() ) { 
    auto counter = 0;
    auto id = G4UIcommand::ConvertToInt(parameters[counter++]); 
    auto title = parameters[counter++];
    fManager->SetH2Title(id, title);     
  }
  else if ( command == fSetH2XAxisCmd.get() ) { 
    auto counter = 0;
    auto id = G4UIcommand::ConvertToInt(parameters[counter++]); 
    auto xaxis = parameters[counter++];
    fManager->SetH2XAxisTitle(id, xaxis);     
  }
  else if ( command == fSetH2YAxisCmd.get() ) { 
    auto counter = 0;
    auto id = G4UIcommand::ConvertToInt(parameters[counter++]); 
    auto yaxis = parameters[counter++];
    fManager->SetH2YAxisTitle(id, yaxis);     
  }
  else if ( command == fSetH2ZAxisCmd.get() ) { 
    auto counter = 0;
    auto id = G4UIcommand::ConvertToInt(parameters[counter++]); 
    auto zaxis = parameters[counter++];
    fManager->SetH2ZAxisTitle(id, zaxis);     
  }
}  
