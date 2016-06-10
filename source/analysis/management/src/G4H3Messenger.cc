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

#include "G4H3Messenger.hh"
#include "G4VAnalysisManager.hh"
#include "G4AnalysisUtilities.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"

#include <iostream>

using namespace G4Analysis;

//_____________________________________________________________________________
G4H3Messenger::G4H3Messenger(G4VAnalysisManager* manager)
  : G4UImessenger(),
    fManager(manager),
    fHelper(nullptr),
    fDirectory(nullptr),
    fCreateH3Cmd(nullptr),
    fSetH3Cmd(nullptr),
    fSetH3XCmd(nullptr),
    fSetH3YCmd(nullptr),
    fSetH3ZCmd(nullptr),
    fSetH3TitleCmd(nullptr), 
    fSetH3XAxisCmd(nullptr), 
    fSetH3YAxisCmd(nullptr), 
    fSetH3ZAxisCmd(nullptr), 
    fXId(-1),
    fYId(-1),
    fXData(),
    fYData()
{  
  fHelper = G4Analysis::make_unique<G4AnalysisMessengerHelper>("h3");

  fDirectory = fHelper->CreateHnDirectory();

  CreateH3Cmd();

  SetH3Cmd();
  fSetH3XCmd = fHelper->CreateSetBinsCommand("x", this);
  fSetH3YCmd = fHelper->CreateSetBinsCommand("y", this);
  
  fSetH3TitleCmd = fHelper->CreateSetTitleCommand(this);
  fSetH3XAxisCmd = fHelper->CreateSetAxisCommand("x", this);
  fSetH3YAxisCmd = fHelper->CreateSetAxisCommand("y", this);
  fSetH3ZAxisCmd = fHelper->CreateSetAxisCommand("z", this);
}

//_____________________________________________________________________________
G4H3Messenger::~G4H3Messenger()
{}

//
// private functions
//

//_____________________________________________________________________________
void G4H3Messenger::CreateH3Cmd()
{
  auto h3Name = new G4UIparameter("name", 's', false);
  h3Name->SetGuidance("Histogram name (label)");
  
  auto h3Title = new G4UIparameter("title", 's', false);
  h3Title->SetGuidance("Histogram title");

  auto h3xNbins0 = new G4UIparameter("xnbins0", 'i', true);
  h3xNbins0->SetGuidance("Number of x-bins (default = 100)");
  h3xNbins0->SetGuidance("Can be reset with /analysis/h3/set command");
  h3xNbins0->SetDefaultValue(100);
  
  auto h3xValMin0 = new G4UIparameter("xvalMin0", 'd', true);
  h3xValMin0->SetGuidance("Minimum x-value, expressed in unit (default = 0.)");
  h3xValMin0->SetGuidance("Can be reset with /analysis/h3/set command");
  h3xValMin0->SetDefaultValue(0.);
  
  auto h3xValMax0 = new G4UIparameter("xvalMax0", 'd', true);
  h3xValMax0->SetGuidance("Maximum x-value, expressed in unit (default = 1.)");
  h3xValMax0->SetGuidance("Can be reset with /analysis/h3/set command");
  h3xValMax0->SetDefaultValue(1.);

  auto h3xValUnit0 = new G4UIparameter("xvalUnit0", 's', true);
  h3xValUnit0->SetGuidance("The unit applied to filled x-values and xvalMin0, xvalMax0");
  h3xValUnit0->SetDefaultValue("none");
  
  auto h3xValFcn0 = new G4UIparameter("xvalFcn0", 's', true);
  G4String fcnxGuidance = "The function applied to filled x-values (log, log10, exp, none).";
  h3xValFcn0->SetGuidance(fcnxGuidance);
  h3xValFcn0->SetParameterCandidates("log log10 exp none");
  h3xValFcn0->SetDefaultValue("none");
  
  auto h3xValBinScheme0 = new G4UIparameter("xvalBinScheme0", 's', true);
  G4String xbinSchemeGuidance = "The binning scheme (linear, log).";
  h3xValBinScheme0->SetParameterCandidates("linear log");
  h3xValBinScheme0->SetGuidance(xbinSchemeGuidance);
  h3xValBinScheme0->SetDefaultValue("linear");
  
  auto h3yNbins0 = new G4UIparameter("ynbins0", 'i', true);
  h3yNbins0->SetGuidance("Number of y-bins (default = 100)");
  h3yNbins0->SetGuidance("Can be reset with /analysis/h3/set command");
  h3yNbins0->SetDefaultValue(100);
  
  auto h3yValMin0 = new G4UIparameter("yvalMin0", 'd', true);
  h3yValMin0->SetGuidance("Minimum y-value, expressed in unit (default = 0.)");
  h3yValMin0->SetGuidance("Can be reset with /analysis/h3/set command");
  h3yValMin0->SetDefaultValue(0.);
  
  auto h3yValMax0 = new G4UIparameter("yvalMax0", 'd', true);
  h3yValMax0->SetGuidance("Maximum y-value, expressed in unit (default = 1.)");
  h3yValMax0->SetGuidance("Can be reset with /analysis/h3/set command");
  h3yValMax0->SetDefaultValue(1.);

  auto h3yValUnit0 = new G4UIparameter("yvalUnit0", 's', true);
  h3yValUnit0->SetGuidance("The unit applied to filled y-values and yvalMin0, yvalMax0");
  h3yValUnit0->SetDefaultValue("none");
  
  auto h3yValFcn0 = new G4UIparameter("yvalFcn0", 's', true);
  G4String fcnyGuidance = "The function applied to filled y-values (log, log10, exp, none).";
  h3yValFcn0->SetGuidance(fcnyGuidance);
  h3yValFcn0->SetParameterCandidates("log log10 exp none");
  h3yValFcn0->SetDefaultValue("none");

  auto h3yValBinScheme0 = new G4UIparameter("yvalBinScheme0", 's', true);
  G4String ybinSchemeGuidance = "The binning scheme (linear, log).";
  h3yValBinScheme0->SetParameterCandidates("linear log");
  h3yValBinScheme0->SetGuidance(ybinSchemeGuidance);
  h3yValBinScheme0->SetDefaultValue("linear");
  
  auto h3zNbins0 = new G4UIparameter("znbins0", 'i', true);
  h3zNbins0->SetGuidance("Number of z-bins (default = 100)");
  h3zNbins0->SetGuidance("Can be reset with /analysis/h3/set command");
  h3zNbins0->SetDefaultValue(100);
  
  auto h3zValMin0 = new G4UIparameter("zvalMin0", 'd', true);
  h3zValMin0->SetGuidance("Minimum z-value, expressed in unit (default = 0.)");
  h3zValMin0->SetGuidance("Can be reset with /analysis/h3/set command");
  h3zValMin0->SetDefaultValue(0.);
  
  auto h3zValMax0 = new G4UIparameter("zvalMax0", 'd', true);
  h3zValMax0->SetGuidance("Maximum z-value, expressed in unit (default = 1.)");
  h3zValMax0->SetGuidance("Can be reset with /analysis/h3/set command");
  h3zValMax0->SetDefaultValue(1.);

  auto h3zValUnit0 = new G4UIparameter("zvalUnit0", 's', true);
  h3zValUnit0->SetGuidance("The unit applied to filled z-values and zvalMin0, zvalMax0");
  h3zValUnit0->SetDefaultValue("none");
  
  auto h3zValFcn0 = new G4UIparameter("zvalFcn0", 's', true);
  G4String fcnzGuidance = "The function applied to filled z-values (log, log10, exp, none).";
  h3zValFcn0->SetGuidance(fcnzGuidance);
  h3zValFcn0->SetParameterCandidates("log log10 exp none");
  h3zValFcn0->SetDefaultValue("none");

  auto h3zValBinScheme0 = new G4UIparameter("zvalBinScheme0", 's', true);
  G4String zbinSchemeGuidance = "The binning scheme (linear, log).";
  h3zValBinScheme0->SetParameterCandidates("linear log");
  h3zValBinScheme0->SetGuidance(zbinSchemeGuidance);
  h3zValBinScheme0->SetDefaultValue("linear");
  
  fCreateH3Cmd = G4Analysis::make_unique<G4UIcommand>("/analysis/h3/create", this);
  fCreateH3Cmd->SetGuidance("Create 3D histogram");
  fCreateH3Cmd->SetParameter(h3Name);
  fCreateH3Cmd->SetParameter(h3Title);
  fCreateH3Cmd->SetParameter(h3xNbins0);
  fCreateH3Cmd->SetParameter(h3xValMin0);
  fCreateH3Cmd->SetParameter(h3xValMax0);
  fCreateH3Cmd->SetParameter(h3xValUnit0);
  fCreateH3Cmd->SetParameter(h3xValFcn0);
  fCreateH3Cmd->SetParameter(h3xValBinScheme0);
  fCreateH3Cmd->SetParameter(h3yNbins0);
  fCreateH3Cmd->SetParameter(h3yValMin0);
  fCreateH3Cmd->SetParameter(h3yValMax0);
  fCreateH3Cmd->SetParameter(h3yValUnit0);
  fCreateH3Cmd->SetParameter(h3yValFcn0);
  fCreateH3Cmd->SetParameter(h3yValBinScheme0);
  fCreateH3Cmd->SetParameter(h3zNbins0);
  fCreateH3Cmd->SetParameter(h3zValMin0);
  fCreateH3Cmd->SetParameter(h3zValMax0);
  fCreateH3Cmd->SetParameter(h3zValUnit0);
  fCreateH3Cmd->SetParameter(h3zValFcn0);
  fCreateH3Cmd->SetParameter(h3zValBinScheme0);
  fCreateH3Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  


//_____________________________________________________________________________
void G4H3Messenger::SetH3Cmd()
{
  auto h3Id = new G4UIparameter("id", 'i', false);
  h3Id->SetGuidance("Histogram id");
  h3Id->SetParameterRange("id>=0");
  
  auto h3xNbins = new G4UIparameter("xnbins", 'i', false);
  h3xNbins->SetGuidance("Number of x-bins");
  
  auto h3xValMin = new G4UIparameter("xvalMin", 'd', false);
  h3xValMin->SetGuidance("Minimum x-value, expressed in unit");
  
  auto h3xValMax = new G4UIparameter("xvalMax", 'd', false);
  h3xValMax->SetGuidance("Maximum x-value, expressed in unit");
  
  auto h3xValUnit = new G4UIparameter("xvalUnit", 's', false);
  h3xValUnit->SetGuidance("The unit applied to filled x-values and xvalMin0, xvalMax0");
  h3xValUnit->SetDefaultValue("none");
 
  auto h3xValFcn = new G4UIparameter("xvalFcn", 's', false);
  h3xValFcn->SetParameterCandidates("log log10 exp none");
  G4String fcnxGuidance = "The function applied to filled x-values (log, log10, exp, none).";
  h3xValFcn->SetGuidance(fcnxGuidance);
  h3xValFcn->SetDefaultValue("none");
 
  auto h3xValBinScheme = new G4UIparameter("xvalBinScheme", 's', true);
  G4String xbinSchemeGuidance = "The binning scheme (linear, log).";
  h3xValBinScheme->SetParameterCandidates("linear log");
  h3xValBinScheme->SetGuidance(xbinSchemeGuidance);
  h3xValBinScheme->SetDefaultValue("linear");
  
  auto h3yNbins = new G4UIparameter("nybins", 'i', false);
  h3yNbins->SetGuidance("Number of y-bins");
  
  auto h3yValMin = new G4UIparameter("yvalMin", 'd', false);
  h3yValMin->SetGuidance("Minimum y-value, expressed in unit");
  
  auto h3yValMax = new G4UIparameter("yvalMax", 'd', false);
  h3yValMax->SetGuidance("Maximum y-value, expressed in unit");
  
  auto h3yValUnit = new G4UIparameter("yvalUnit", 's', true);
  h3yValUnit->SetGuidance("The unit applied to filled y-values and yvalMin0, yvalMax0");
  h3yValUnit->SetDefaultValue("none");
 
  auto h3yValFcn = new G4UIparameter("yvalFcn", 's', false);
  h3yValFcn->SetParameterCandidates("log log10 exp none");
  G4String fcnyGuidance = "The function applied to filled y-values (log, log10, exp, none).";
  h3yValFcn->SetGuidance(fcnyGuidance);
  h3yValFcn->SetDefaultValue("none");
 
  auto h3yValBinScheme = new G4UIparameter("yvalBinScheme", 's', true);
  G4String ybinSchemeGuidance = "The binning scheme (linear, log).";
  h3yValBinScheme->SetParameterCandidates("linear log");
  h3yValBinScheme->SetGuidance(ybinSchemeGuidance);
  h3yValBinScheme->SetDefaultValue("linear");
 
  auto h3zNbins = new G4UIparameter("nzbins", 'i', false);
  h3zNbins->SetGuidance("Number of z-bins");
  
  auto h3zValMin = new G4UIparameter("zvalMin", 'd', false);
  h3zValMin->SetGuidance("Minimum z-value, expressed in unit");
  
  auto h3zValMax = new G4UIparameter("zvalMax", 'd', false);
  h3zValMax->SetGuidance("Maximum z-value, expressed in unit");
  
  auto h3zValUnit = new G4UIparameter("zvalUnit", 's', true);
  h3zValUnit->SetGuidance("The unit applied to filled z-values and zvalMin0, zvalMax0");
  h3zValUnit->SetDefaultValue("none");
 
  auto h3zValFcn = new G4UIparameter("zvalFcn", 's', false);
  h3zValFcn->SetParameterCandidates("log log10 exp none");
  G4String fcnzGuidance = "The function applied to filled z-values (log, log10, exp, none).";
  h3zValFcn->SetGuidance(fcnzGuidance);
  h3zValFcn->SetDefaultValue("none");
 
  auto h3zValBinScheme = new G4UIparameter("zvalBinScheme", 's', true);
  G4String zbinSchemeGuidance = "The binning scheme (linear, log).";
  h3zValBinScheme->SetParameterCandidates("linear log");
  h3zValBinScheme->SetGuidance(zbinSchemeGuidance);
  h3zValBinScheme->SetDefaultValue("linear");
 
  fSetH3Cmd = G4Analysis::make_unique<G4UIcommand>("/analysis/h3/set", this);
  fSetH3Cmd->SetGuidance("Set parameters for the 3D histogram of given id:");
  fSetH3Cmd->SetGuidance("  nxbins; xvalMin; xvalMax; xunit; xfunction; xbinScheme");
  fSetH3Cmd->SetGuidance("  nybins; yvalMin; yvalMax; yunit; yfunction; ybinScheme");
  fSetH3Cmd->SetGuidance("  nzbins; zvalMin; zvalMax; zunit; zfunction; zbinScheme");
  fSetH3Cmd->SetParameter(h3Id);
  fSetH3Cmd->SetParameter(h3xNbins);
  fSetH3Cmd->SetParameter(h3xValMin);
  fSetH3Cmd->SetParameter(h3xValMax);
  fSetH3Cmd->SetParameter(h3xValUnit);
  fSetH3Cmd->SetParameter(h3xValFcn);
  fSetH3Cmd->SetParameter(h3xValBinScheme);
  fSetH3Cmd->SetParameter(h3yNbins);
  fSetH3Cmd->SetParameter(h3yValMin);
  fSetH3Cmd->SetParameter(h3yValMax);
  fSetH3Cmd->SetParameter(h3yValUnit);
  fSetH3Cmd->SetParameter(h3yValFcn);
  fSetH3Cmd->SetParameter(h3yValBinScheme);
  fSetH3Cmd->SetParameter(h3zNbins);
  fSetH3Cmd->SetParameter(h3zValMin);
  fSetH3Cmd->SetParameter(h3zValMax);
  fSetH3Cmd->SetParameter(h3zValUnit);
  fSetH3Cmd->SetParameter(h3zValFcn);
  fSetH3Cmd->SetParameter(h3zValBinScheme);
  fSetH3Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  


//
// public functions
//

//_____________________________________________________________________________
void G4H3Messenger::SetNewValue(G4UIcommand* command, G4String newValues)
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

  if ( command == fCreateH3Cmd.get() ) {  
    auto counter = 0;
    auto name = parameters[counter++];
    auto title = parameters[counter++];
    G4AnalysisMessengerHelper::BinData xdata;
    fHelper->GetBinData(xdata, parameters, counter);
    auto xunit = GetUnitValue(xdata.fSunit);
    G4AnalysisMessengerHelper::BinData ydata;
    fHelper->GetBinData(ydata, parameters, counter);
    auto yunit = GetUnitValue(ydata.fSunit);
    G4AnalysisMessengerHelper::BinData zdata;
    fHelper->GetBinData(zdata, parameters, counter);
    auto zunit = GetUnitValue(zdata.fSunit);
    fManager->CreateH3(name, title, 
                       xdata.fNbins, xdata.fVmin*xunit, xdata.fVmax*xunit, 
                       ydata.fNbins, ydata.fVmin*yunit, ydata.fVmax*yunit, 
                       zdata.fNbins, zdata.fVmin*zunit, zdata.fVmax*zunit, 
                       xdata.fSunit, ydata.fSunit, zdata.fSunit,
                       xdata.fSfcn, ydata.fSfcn, zdata.fSfcn,
                       xdata.fSbinScheme, ydata.fSbinScheme, zdata.fSbinScheme);     
  }
  else if ( command == fSetH3Cmd.get() ) { 
    auto counter = 0;
    auto id = G4UIcommand::ConvertToInt(parameters[counter++]);
    G4AnalysisMessengerHelper::BinData xdata;
    fHelper->GetBinData(xdata, parameters, counter);
    auto xunit = GetUnitValue(xdata.fSunit);
    G4AnalysisMessengerHelper::BinData ydata;
    fHelper->GetBinData(ydata, parameters, counter);
    auto yunit = GetUnitValue(ydata.fSunit);
    G4AnalysisMessengerHelper::BinData zdata;
    fHelper->GetBinData(zdata, parameters, counter);
    auto zunit = GetUnitValue(zdata.fSunit);
    fManager->SetH3(id, 
                    xdata.fNbins, xdata.fVmin*xunit, xdata.fVmax*xunit, 
                    ydata.fNbins, ydata.fVmin*yunit, ydata.fVmax*yunit, 
                    zdata.fNbins, zdata.fVmin*zunit, zdata.fVmax*zunit, 
                    xdata.fSunit, ydata.fSunit, zdata.fSunit,
                    xdata.fSfcn, ydata.fSfcn, zdata.fSfcn,
                    xdata.fSbinScheme, ydata.fSbinScheme, zdata.fSbinScheme);     
  }
  else if ( command == fSetH3XCmd.get() ) { 
    // Only save values
    auto counter = 0;
    fXId = G4UIcommand::ConvertToInt(parameters[counter++]);
    fHelper->GetBinData(fXData, parameters, counter);
  }
  else if ( command == fSetH3YCmd.get() ) { 
    // Only save values
    auto counter = 0;
    fYId = G4UIcommand::ConvertToInt(parameters[counter++]);
    fHelper->GetBinData(fYData, parameters, counter);
  }
  else if ( command == fSetH3ZCmd.get() ) { 
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
    G4AnalysisMessengerHelper::BinData zdata;
    fHelper->GetBinData(zdata, parameters, counter);
    auto zunit = GetUnitValue(zdata.fSunit);
    fManager->SetH3(id, 
                    fXData.fNbins, fXData.fVmin*xunit, fXData.fVmax*xunit,
                    fYData.fNbins, fYData.fVmin*yunit, fYData.fVmax*yunit,
                    zdata.fNbins, zdata.fVmin*zunit, zdata.fVmax*zunit, 
                    fXData.fSunit, fYData.fSunit, zdata.fSunit, 
                    fXData.fSfcn, fYData.fSfcn, zdata.fSfcn,
                    fXData.fSbinScheme, fYData.fSbinScheme, zdata.fSbinScheme);     
    fXId = -1;                        
    fYId = -1;                        
  }
  else if ( command == fSetH3TitleCmd.get() ) { 
    auto counter = 0;
    auto id = G4UIcommand::ConvertToInt(parameters[counter++]); 
    auto title = parameters[counter++];
    fManager->SetH3Title(id, title);     
  }
  else if ( command == fSetH3XAxisCmd.get() ) { 
    auto counter = 0;
    auto id = G4UIcommand::ConvertToInt(parameters[counter++]); 
    auto xaxis = parameters[counter++];
    fManager->SetH3XAxisTitle(id, xaxis);     
  }
  else if ( command == fSetH3YAxisCmd.get() ) { 
    auto counter = 0;
    auto id = G4UIcommand::ConvertToInt(parameters[counter++]); 
    auto yaxis = parameters[counter++];
    fManager->SetH3YAxisTitle(id, yaxis);     
  }
  else if ( command == fSetH3ZAxisCmd.get() ) { 
    auto counter = 0;
    auto id = G4UIcommand::ConvertToInt(parameters[counter++]); 
    auto zaxis = parameters[counter++];
    fManager->SetH3ZAxisTitle(id, zaxis);     
  }
}  
