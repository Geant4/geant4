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

namespace {

void Exception(G4UIcommand* command, G4int nofParameters)
{
  G4ExceptionDescription description;
  description 
    << "Got wrong number of \"" << command->GetCommandName() 
    << "\" parameters: " << nofParameters
    << " instead of " << command->GetParameterEntries() 
    << " expected" << G4endl;
  G4Exception("G4H2Messenger::SetNewValue",
              "Analysis_W013", JustWarning, description);
}

}                  


//_____________________________________________________________________________
G4H2Messenger::G4H2Messenger(G4VAnalysisManager* manager)
  : G4UImessenger(),
    fManager(manager),
    fH2Dir(0),  
    fCreateH2Cmd(0),
    fSetH2Cmd(0),
    fSetH2TitleCmd(0), 
    fSetH2XAxisCmd(0), 
    fSetH2YAxisCmd(0)
{  
  fH2Dir = new G4UIdirectory("/analysis/h2/");
  fH2Dir->SetGuidance("2D histograms control");

  CreateH2Cmd();
  SetH2Cmd();
  
  SetH2TitleCmd();
  SetH2XAxisCmd();
  SetH2YAxisCmd();
  SetH2ZAxisCmd();
}

//_____________________________________________________________________________
G4H2Messenger::~G4H2Messenger()
{
  delete fCreateH2Cmd;
  delete fSetH2Cmd;
  delete fSetH2TitleCmd;  
  delete fSetH2XAxisCmd;  
  delete fSetH2YAxisCmd;  
  delete fSetH2ZAxisCmd;  
  delete fH2Dir;
}

//
// private functions
//

//_____________________________________________________________________________
void G4H2Messenger::CreateH2Cmd()
{
  G4UIparameter* h2Name = new G4UIparameter("name", 's', false);
  h2Name->SetGuidance("Histogram name (label)");
  
  G4UIparameter* h2Title = new G4UIparameter("title", 's', false);
  h2Title->SetGuidance("Histogram title");

  G4UIparameter* h2xNbins0 = new G4UIparameter("xnbins0", 'i', true);
  h2xNbins0->SetGuidance("Number of x-bins (default = 100)");
  h2xNbins0->SetGuidance("Can be reset with /analysis/h2/set command");
  h2xNbins0->SetDefaultValue(100);
  
  G4UIparameter* h2xValMin0 = new G4UIparameter("xvalMin0", 'd', true);
  h2xValMin0->SetGuidance("Minimum x-value, expressed in unit (default = 0.)");
  h2xValMin0->SetGuidance("Can be reset with /analysis/h2/set command");
  h2xValMin0->SetDefaultValue(0.);
  
  G4UIparameter* h2xValMax0 = new G4UIparameter("xvalMax0", 'd', true);
  h2xValMax0->SetGuidance("Maximum x-value, expressed in unit (default = 1.)");
  h2xValMax0->SetGuidance("Can be reset with /analysis/h2/set command");
  h2xValMax0->SetDefaultValue(1.);

  G4UIparameter* h2xValUnit0 = new G4UIparameter("xvalUnit0", 's', true);
  h2xValUnit0->SetGuidance("The unit applied to filled x-values and xvalMin0, xvalMax0");
  h2xValUnit0->SetDefaultValue("none");
  
  G4UIparameter* h2xValFcn0 = new G4UIparameter("xvalFcn0", 's', true);
  G4String fcnxGuidance = "The function applied to filled x-values (log, log10, exp, none).";
  h2xValFcn0->SetGuidance(fcnxGuidance);
  h2xValFcn0->SetParameterCandidates("log log10 exp none");
  h2xValFcn0->SetDefaultValue("none");
  
  G4UIparameter* h2xValBinScheme0 = new G4UIparameter("xvalBinScheme0", 's', true);
  G4String xbinSchemeGuidance = "The binning scheme (linear, log).";
  h2xValBinScheme0->SetParameterCandidates("linear log");
  h2xValBinScheme0->SetGuidance(xbinSchemeGuidance);
  h2xValBinScheme0->SetDefaultValue("linear");
  
  G4UIparameter* h2yNbins0 = new G4UIparameter("ynbins0", 'i', true);
  h2yNbins0->SetGuidance("Number of y-bins (default = 100)");
  h2yNbins0->SetGuidance("Can be reset with /analysis/h2/set command");
  h2yNbins0->SetDefaultValue(100);
  
  G4UIparameter* h2yValMin0 = new G4UIparameter("yvalMin0", 'd', true);
  h2yValMin0->SetGuidance("Minimum y-value, expressed in unit (default = 0.)");
  h2yValMin0->SetGuidance("Can be reset with /analysis/h2/set command");
  h2yValMin0->SetDefaultValue(0.);
  
  G4UIparameter* h2yValMax0 = new G4UIparameter("yvalMax0", 'd', true);
  h2yValMax0->SetGuidance("Maximum y-value, expressed in unit (default = 1.)");
  h2yValMax0->SetGuidance("Can be reset with /analysis/h2/set command");
  h2yValMax0->SetDefaultValue(1.);

  G4UIparameter* h2yValUnit0 = new G4UIparameter("yvalUnit0", 's', true);
  h2yValUnit0->SetGuidance("The unit applied to filled y-values and yvalMin0, yvalMax0");
  h2yValUnit0->SetDefaultValue("none");
  
  G4UIparameter* h2yValFcn0 = new G4UIparameter("yvalFcn0", 's', true);
  G4String fcnyGuidance = "The function applied to filled y-values (log, log10, exp, none).";
  h2yValFcn0->SetGuidance(fcnyGuidance);
  h2yValFcn0->SetParameterCandidates("log log10 exp none");
  h2yValFcn0->SetDefaultValue("none");

  G4UIparameter* h2yValBinScheme0 = new G4UIparameter("yvalBinScheme0", 's', true);
  G4String ybinSchemeGuidance = "The binning scheme (linear, log).";
  h2yValBinScheme0->SetParameterCandidates("linear log");
  h2yValBinScheme0->SetGuidance(ybinSchemeGuidance);
  h2yValBinScheme0->SetDefaultValue("linear");
  
  fCreateH2Cmd = new G4UIcommand("/analysis/h2/create", this);
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
  G4UIparameter* h2Id = new G4UIparameter("id", 'i', false);
  h2Id->SetGuidance("Histogram id");
  h2Id->SetParameterRange("id>=0");
  
  G4UIparameter* h2xNbins = new G4UIparameter("xnbins", 'i', false);
  h2xNbins->SetGuidance("Number of x-bins");
  
  G4UIparameter* h2xValMin = new G4UIparameter("xvalMin", 'd', false);
  h2xValMin->SetGuidance("Minimum x-value, expressed in unit");
  
  G4UIparameter* h2xValMax = new G4UIparameter("xvalMax", 'd', false);
  h2xValMax->SetGuidance("Maximum x-value, expressed in unit");
  
  G4UIparameter* h2xValUnit = new G4UIparameter("xvalUnit", 's', false);
  h2xValUnit->SetGuidance("The unit applied to filled x-values and xvalMin, xvalMax");
  h2xValUnit->SetDefaultValue("none");
 
  G4UIparameter* h2xValFcn = new G4UIparameter("xvalFcn", 's', false);
  h2xValFcn->SetParameterCandidates("log log10 exp none");
  G4String fcnxGuidance = "The function applied to filled x-values (log, log10, exp, none).";
  h2xValFcn->SetGuidance(fcnxGuidance);
  h2xValFcn->SetDefaultValue("none");
 
  G4UIparameter* h2xValBinScheme = new G4UIparameter("xvalBinScheme", 's', true);
  G4String xbinSchemeGuidance = "The binning scheme (linear, log).";
  h2xValBinScheme->SetParameterCandidates("linear log");
  h2xValBinScheme->SetGuidance(xbinSchemeGuidance);
  h2xValBinScheme->SetDefaultValue("linear");
  
  G4UIparameter* h2yNbins = new G4UIparameter("nybins", 'i', false);
  h2yNbins->SetGuidance("Number of y-bins");
  
  G4UIparameter* h2yValMin = new G4UIparameter("yvalMin", 'd', false);
  h2yValMin->SetGuidance("Minimum y-value, expressed in unit");
  
  G4UIparameter* h2yValMax = new G4UIparameter("yvalMax", 'd', false);
  h2yValMax->SetGuidance("Maximum y-value, expressed in unit");
  
  G4UIparameter* h2yValUnit = new G4UIparameter("yvalUnit", 's', true);
  h2yValUnit->SetGuidance("The unit applied to filled y-values and yvalMin, yvalMax");
  h2yValUnit->SetDefaultValue("none");
 
  G4UIparameter* h2yValFcn = new G4UIparameter("yvalFcn", 's', false);
  h2yValFcn->SetParameterCandidates("log log10 exp none");
  G4String fcnyGuidance = "The function applied to filled y-values (log, log10, exp, none).";
  h2yValFcn->SetGuidance(fcnyGuidance);
  h2yValFcn->SetDefaultValue("none");
 
  G4UIparameter* h2yValBinScheme = new G4UIparameter("yvalBinScheme", 's', true);
  G4String ybinSchemeGuidance = "The binning scheme (linear, log).";
  h2yValBinScheme->SetParameterCandidates("linear log");
  h2yValBinScheme->SetGuidance(ybinSchemeGuidance);
  h2yValBinScheme->SetDefaultValue("linear");

  fSetH2Cmd = new G4UIcommand("/analysis/h2/set", this);
  fSetH2Cmd->SetGuidance("Set parameters for the 2D histogram of #Id :");
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

//_____________________________________________________________________________
void G4H2Messenger::SetH2TitleCmd()
{
  G4UIparameter* h2Id = new G4UIparameter("idTitle", 'i', false);
  h2Id->SetGuidance("Histogram id");
  h2Id->SetParameterRange("idTitle>=0");

  G4UIparameter* h2Title = new G4UIparameter("h2Title", 's', true);
  h2Title->SetGuidance("Histogram title");
  h2Title->SetDefaultValue("none");

  fSetH2TitleCmd = new G4UIcommand("/analysis/h2/setTitle", this);
  fSetH2TitleCmd->SetGuidance("Set title for the 2D histogram of #Id");
  fSetH2TitleCmd->SetParameter(h2Id);
  fSetH2TitleCmd->SetParameter(h2Title);
  fSetH2TitleCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  

//_____________________________________________________________________________
void G4H2Messenger::SetH2XAxisCmd()
{
  G4UIparameter* h2Id = new G4UIparameter("idXaxis", 'i', false);
  h2Id->SetGuidance("Histogram id");
  h2Id->SetParameterRange("idXaxis>=0");

  G4UIparameter* h2XAxis = new G4UIparameter("h2Xaxis", 's', true);
  h2XAxis->SetGuidance("Histogram x-axis title");
  h2XAxis->SetDefaultValue("none");

  fSetH2XAxisCmd = new G4UIcommand("/analysis/h2/setXaxis", this);
  fSetH2XAxisCmd->SetGuidance("Set x-axis title for the 2D histogram of #Id");
  fSetH2XAxisCmd->SetParameter(h2Id);
  fSetH2XAxisCmd->SetParameter(h2XAxis);
  fSetH2XAxisCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  

//_____________________________________________________________________________
void G4H2Messenger::SetH2YAxisCmd()
{
  G4UIparameter* h2Id = new G4UIparameter("idYaxis", 'i', false);
  h2Id->SetGuidance("Histogram id");
  h2Id->SetParameterRange("idYaxis>=0");

  G4UIparameter* h2YAxis = new G4UIparameter("h2Yaxis", 's', true);
  h2YAxis->SetGuidance("Histogram y-axis title");
  h2YAxis->SetDefaultValue("none");

  fSetH2YAxisCmd = new G4UIcommand("/analysis/h2/setYaxis", this);
  fSetH2YAxisCmd->SetGuidance("Set y-axis title for the 2D histogram of #Id");
  fSetH2YAxisCmd->SetParameter(h2Id);
  fSetH2YAxisCmd->SetParameter(h2YAxis);
  fSetH2YAxisCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  

//_____________________________________________________________________________
void G4H2Messenger::SetH2ZAxisCmd()
{
  G4UIparameter* h2Id = new G4UIparameter("idYaxis", 'i', false);
  h2Id->SetGuidance("Histogram id");
  h2Id->SetParameterRange("idYaxis>=0");

  G4UIparameter* h2YAxis = new G4UIparameter("h2Yaxis", 's', true);
  h2YAxis->SetGuidance("Histogram y-axis title");
  h2YAxis->SetDefaultValue("none");

  fSetH2ZAxisCmd = new G4UIcommand("/analysis/h2/setYaxis", this);
  fSetH2ZAxisCmd->SetGuidance("Set y-axis title for the 2D histogram of #Id");
  fSetH2ZAxisCmd->SetParameter(h2Id);
  fSetH2ZAxisCmd->SetParameter(h2YAxis);
  fSetH2ZAxisCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
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
    Exception(command, parameters.size());
    return;
  }  

  if ( command == fCreateH2Cmd ) { 
    G4int counter = 0;
    G4String name = parameters[counter++];
    G4String title = parameters[counter++];
    G4int xnbins = G4UIcommand::ConvertToInt(parameters[counter++]); 
    G4double xvmin = G4UIcommand::ConvertToDouble(parameters[counter++]); 
    G4double xvmax = G4UIcommand::ConvertToDouble(parameters[counter++]); ; 
    G4String xsunit = parameters[counter++];
    G4String xsfcn = parameters[counter++];
    G4String xsbinScheme = parameters[counter++];
    G4double xunit = GetUnitValue(xsunit);
    G4int ynbins = G4UIcommand::ConvertToInt(parameters[counter++]); 
    G4double yvmin = G4UIcommand::ConvertToDouble(parameters[counter++]); 
    G4double yvmax = G4UIcommand::ConvertToDouble(parameters[counter++]); ; 
    G4String ysunit = parameters[counter++];
    G4String ysfcn = parameters[counter++];
    G4String ysbinScheme = parameters[counter++];
    G4double yunit = GetUnitValue(ysunit);
    fManager->CreateH2(name, title, 
                       xnbins, xvmin*xunit, xvmax*xunit, 
                       ynbins, yvmin*yunit, yvmax*yunit, 
                       xsunit, ysunit, xsfcn, ysfcn, xsbinScheme, ysbinScheme);     
  }
  else if ( command == fSetH2Cmd ) {
    G4int counter = 0;
    G4int id = G4UIcommand::ConvertToInt(parameters[counter++]);
    G4int xnbins = G4UIcommand::ConvertToInt(parameters[counter++]); 
    G4double xvmin = G4UIcommand::ConvertToDouble(parameters[counter++]); 
    G4double xvmax = G4UIcommand::ConvertToDouble(parameters[counter++]); ; 
    G4String xsunit = parameters[counter++];
    G4String xsfcn = parameters[counter++];
    G4String xsbinScheme = parameters[counter++];
    G4double xunit = GetUnitValue(xsunit);
    G4int ynbins = G4UIcommand::ConvertToInt(parameters[counter++]); 
    G4double yvmin = G4UIcommand::ConvertToDouble(parameters[counter++]); 
    G4double yvmax = G4UIcommand::ConvertToDouble(parameters[counter++]); ; 
    G4String ysunit = parameters[counter++];
    G4String ysfcn = parameters[counter++];
    G4String ysbinScheme = parameters[counter++];
    G4double yunit = GetUnitValue(ysunit);
    fManager->SetH2(id, 
                    xnbins, xvmin*xunit, xvmax*xunit,
                    ynbins, yvmin*yunit, yvmax*yunit, 
                    xsunit, ysunit, xsfcn, ysfcn, xsbinScheme, ysbinScheme);     
  }
  else if ( command == fSetH2TitleCmd ) {
    G4int counter = 0;
    G4int id = G4UIcommand::ConvertToInt(parameters[counter++]); 
    G4String title = parameters[counter++];
    fManager->SetH2Title(id, title);     
  }
  else if ( command == fSetH2XAxisCmd ) {
    G4int counter = 0;
    G4int id = G4UIcommand::ConvertToInt(parameters[counter++]); 
    G4String xaxis = parameters[counter++];
    fManager->SetH2XAxisTitle(id, xaxis);     
  }
  else if ( command == fSetH2YAxisCmd ) {
    G4int counter = 0;
    G4int id = G4UIcommand::ConvertToInt(parameters[counter++]); 
    G4String yaxis = parameters[counter++];
    fManager->SetH2YAxisTitle(id, yaxis);     
  }
  else if ( command == fSetH2ZAxisCmd ) {
    G4int counter = 0;
    G4int id = G4UIcommand::ConvertToInt(parameters[counter++]); 
    G4String zaxis = parameters[counter++];
    fManager->SetH2ZAxisTitle(id, zaxis);     
  }
}  
