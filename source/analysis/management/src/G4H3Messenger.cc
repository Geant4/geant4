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

namespace {

void Exception(G4UIcommand* command, G4int nofParameters)
{
  G4ExceptionDescription description;
  description 
    << "Got wrong number of \"" << command->GetCommandName() 
    << "\" parameters: " << nofParameters
    << " instead of " << command->GetParameterEntries() 
    << " expected" << G4endl;
  G4Exception("G4H3Messenger::SetNewValue",
              "Analysis_W013", JustWarning, description);
}

}                  


//_____________________________________________________________________________
G4H3Messenger::G4H3Messenger(G4VAnalysisManager* manager)
  : G4UImessenger(),
    fManager(manager),
    fH3Dir(0),  
    fCreateH3Cmd(0),
    fSetH3Cmd(0),
    fSetH3TitleCmd(0), 
    fSetH3XAxisCmd(0), 
    fSetH3YAxisCmd(0)
{  
  fH3Dir = new G4UIdirectory("/analysis/h3/");
  fH3Dir->SetGuidance("3D histograms control");

  CreateH3Cmd();
  SetH3Cmd();
  
  SetH3TitleCmd();
  SetH3XAxisCmd();
  SetH3YAxisCmd();
  SetH3ZAxisCmd();
}

//_____________________________________________________________________________
G4H3Messenger::~G4H3Messenger()
{
  delete fCreateH3Cmd;
  delete fSetH3Cmd;
  delete fSetH3TitleCmd;  
  delete fSetH3XAxisCmd;  
  delete fSetH3YAxisCmd;  
  delete fSetH3ZAxisCmd;  
  delete fH3Dir;
}

//
// private functions
//

//_____________________________________________________________________________
void G4H3Messenger::CreateH3Cmd()
{
  G4UIparameter* h3Name = new G4UIparameter("name", 's', false);
  h3Name->SetGuidance("Histogram name (label)");
  
  G4UIparameter* h3Title = new G4UIparameter("title", 's', false);
  h3Title->SetGuidance("Histogram title");

  G4UIparameter* h3xNbins0 = new G4UIparameter("xnbins0", 'i', true);
  h3xNbins0->SetGuidance("Number of x-bins (default = 100)");
  h3xNbins0->SetGuidance("Can be reset with /analysis/h3/set command");
  h3xNbins0->SetDefaultValue(100);
  
  G4UIparameter* h3xValMin0 = new G4UIparameter("xvalMin0", 'd', true);
  h3xValMin0->SetGuidance("Minimum x-value, expressed in unit (default = 0.)");
  h3xValMin0->SetGuidance("Can be reset with /analysis/h3/set command");
  h3xValMin0->SetDefaultValue(0.);
  
  G4UIparameter* h3xValMax0 = new G4UIparameter("xvalMax0", 'd', true);
  h3xValMax0->SetGuidance("Maximum x-value, expressed in unit (default = 1.)");
  h3xValMax0->SetGuidance("Can be reset with /analysis/h3/set command");
  h3xValMax0->SetDefaultValue(1.);

  G4UIparameter* h3xValUnit0 = new G4UIparameter("xvalUnit0", 's', true);
  h3xValUnit0->SetGuidance("The unit applied to filled x-values and xvalMin0, xvalMax0");
  h3xValUnit0->SetDefaultValue("none");
  
  G4UIparameter* h3xValFcn0 = new G4UIparameter("xvalFcn0", 's', true);
  G4String fcnxGuidance = "The function applied to filled x-values (log, log10, exp, none).";
  h3xValFcn0->SetGuidance(fcnxGuidance);
  h3xValFcn0->SetParameterCandidates("log log10 exp none");
  h3xValFcn0->SetDefaultValue("none");
  
  G4UIparameter* h3xValBinScheme0 = new G4UIparameter("xvalBinScheme0", 's', true);
  G4String xbinSchemeGuidance = "The binning scheme (linear, log).";
  h3xValBinScheme0->SetParameterCandidates("linear log");
  h3xValBinScheme0->SetGuidance(xbinSchemeGuidance);
  h3xValBinScheme0->SetDefaultValue("linear");
  
  G4UIparameter* h3yNbins0 = new G4UIparameter("ynbins0", 'i', true);
  h3yNbins0->SetGuidance("Number of y-bins (default = 100)");
  h3yNbins0->SetGuidance("Can be reset with /analysis/h3/set command");
  h3yNbins0->SetDefaultValue(100);
  
  G4UIparameter* h3yValMin0 = new G4UIparameter("yvalMin0", 'd', true);
  h3yValMin0->SetGuidance("Minimum y-value, expressed in unit (default = 0.)");
  h3yValMin0->SetGuidance("Can be reset with /analysis/h3/set command");
  h3yValMin0->SetDefaultValue(0.);
  
  G4UIparameter* h3yValMax0 = new G4UIparameter("yvalMax0", 'd', true);
  h3yValMax0->SetGuidance("Maximum y-value, expressed in unit (default = 1.)");
  h3yValMax0->SetGuidance("Can be reset with /analysis/h3/set command");
  h3yValMax0->SetDefaultValue(1.);

  G4UIparameter* h3yValUnit0 = new G4UIparameter("yvalUnit0", 's', true);
  h3yValUnit0->SetGuidance("The unit applied to filled y-values and yvalMin0, yvalMax0");
  h3yValUnit0->SetDefaultValue("none");
  
  G4UIparameter* h3yValFcn0 = new G4UIparameter("yvalFcn0", 's', true);
  G4String fcnyGuidance = "The function applied to filled y-values (log, log10, exp, none).";
  h3yValFcn0->SetGuidance(fcnyGuidance);
  h3yValFcn0->SetParameterCandidates("log log10 exp none");
  h3yValFcn0->SetDefaultValue("none");

  G4UIparameter* h3yValBinScheme0 = new G4UIparameter("yvalBinScheme0", 's', true);
  G4String ybinSchemeGuidance = "The binning scheme (linear, log).";
  h3yValBinScheme0->SetParameterCandidates("linear log");
  h3yValBinScheme0->SetGuidance(ybinSchemeGuidance);
  h3yValBinScheme0->SetDefaultValue("linear");
  
  G4UIparameter* h3zNbins0 = new G4UIparameter("znbins0", 'i', true);
  h3zNbins0->SetGuidance("Number of z-bins (default = 100)");
  h3zNbins0->SetGuidance("Can be reset with /analysis/h3/set command");
  h3zNbins0->SetDefaultValue(100);
  
  G4UIparameter* h3zValMin0 = new G4UIparameter("zvalMin0", 'd', true);
  h3zValMin0->SetGuidance("Minimum z-value, expressed in unit (default = 0.)");
  h3zValMin0->SetGuidance("Can be reset with /analysis/h3/set command");
  h3zValMin0->SetDefaultValue(0.);
  
  G4UIparameter* h3zValMax0 = new G4UIparameter("zvalMax0", 'd', true);
  h3zValMax0->SetGuidance("Maximum z-value, expressed in unit (default = 1.)");
  h3zValMax0->SetGuidance("Can be reset with /analysis/h3/set command");
  h3zValMax0->SetDefaultValue(1.);

  G4UIparameter* h3zValUnit0 = new G4UIparameter("zvalUnit0", 's', true);
  h3zValUnit0->SetGuidance("The unit applied to filled z-values and zvalMin0, zvalMax0");
  h3zValUnit0->SetDefaultValue("none");
  
  G4UIparameter* h3zValFcn0 = new G4UIparameter("zvalFcn0", 's', true);
  G4String fcnzGuidance = "The function applied to filled z-values (log, log10, exp, none).";
  h3zValFcn0->SetGuidance(fcnzGuidance);
  h3zValFcn0->SetParameterCandidates("log log10 exp none");
  h3zValFcn0->SetDefaultValue("none");

  G4UIparameter* h3zValBinScheme0 = new G4UIparameter("zvalBinScheme0", 's', true);
  G4String zbinSchemeGuidance = "The binning scheme (linear, log).";
  h3zValBinScheme0->SetParameterCandidates("linear log");
  h3zValBinScheme0->SetGuidance(zbinSchemeGuidance);
  h3zValBinScheme0->SetDefaultValue("linear");
  
  fCreateH3Cmd = new G4UIcommand("/analysis/h3/create", this);
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
  G4UIparameter* h3Id = new G4UIparameter("id", 'i', false);
  h3Id->SetGuidance("Histogram id");
  h3Id->SetParameterRange("id>=0");
  
  G4UIparameter* h3xNbins = new G4UIparameter("xnbins", 'i', false);
  h3xNbins->SetGuidance("Number of x-bins");
  
  G4UIparameter* h3xValMin = new G4UIparameter("xvalMin", 'd', false);
  h3xValMin->SetGuidance("Minimum x-value, expressed in unit");
  
  G4UIparameter* h3xValMax = new G4UIparameter("xvalMax", 'd', false);
  h3xValMax->SetGuidance("Maximum x-value, expressed in unit");
  
  G4UIparameter* h3xValUnit = new G4UIparameter("xvalUnit", 's', false);
  h3xValUnit->SetGuidance("The unit applied to filled x-values and xvalMin0, xvalMax0");
  h3xValUnit->SetDefaultValue("none");
 
  G4UIparameter* h3xValFcn = new G4UIparameter("xvalFcn", 's', false);
  h3xValFcn->SetParameterCandidates("log log10 exp none");
  G4String fcnxGuidance = "The function applied to filled x-values (log, log10, exp, none).";
  h3xValFcn->SetGuidance(fcnxGuidance);
  h3xValFcn->SetDefaultValue("none");
 
  G4UIparameter* h3xValBinScheme = new G4UIparameter("xvalBinScheme", 's', true);
  G4String xbinSchemeGuidance = "The binning scheme (linear, log).";
  h3xValBinScheme->SetParameterCandidates("linear log");
  h3xValBinScheme->SetGuidance(xbinSchemeGuidance);
  h3xValBinScheme->SetDefaultValue("linear");
  
  G4UIparameter* h3yNbins = new G4UIparameter("nybins", 'i', false);
  h3yNbins->SetGuidance("Number of y-bins");
  
  G4UIparameter* h3yValMin = new G4UIparameter("yvalMin", 'd', false);
  h3yValMin->SetGuidance("Minimum y-value, expressed in unit");
  
  G4UIparameter* h3yValMax = new G4UIparameter("yvalMax", 'd', false);
  h3yValMax->SetGuidance("Maximum y-value, expressed in unit");
  
  G4UIparameter* h3yValUnit = new G4UIparameter("yvalUnit", 's', true);
  h3yValUnit->SetGuidance("The unit applied to filled y-values and yvalMin0, yvalMax0");
  h3yValUnit->SetDefaultValue("none");
 
  G4UIparameter* h3yValFcn = new G4UIparameter("yvalFcn", 's', false);
  h3yValFcn->SetParameterCandidates("log log10 exp none");
  G4String fcnyGuidance = "The function applied to filled y-values (log, log10, exp, none).";
  h3yValFcn->SetGuidance(fcnyGuidance);
  h3yValFcn->SetDefaultValue("none");
 
  G4UIparameter* h3yValBinScheme = new G4UIparameter("yvalBinScheme", 's', true);
  G4String ybinSchemeGuidance = "The binning scheme (linear, log).";
  h3yValBinScheme->SetParameterCandidates("linear log");
  h3yValBinScheme->SetGuidance(ybinSchemeGuidance);
  h3yValBinScheme->SetDefaultValue("linear");
 
  G4UIparameter* h3zNbins = new G4UIparameter("nzbins", 'i', false);
  h3zNbins->SetGuidance("Number of z-bins");
  
  G4UIparameter* h3zValMin = new G4UIparameter("zvalMin", 'd', false);
  h3zValMin->SetGuidance("Minimum z-value, expressed in unit");
  
  G4UIparameter* h3zValMax = new G4UIparameter("zvalMax", 'd', false);
  h3zValMax->SetGuidance("Maximum z-value, expressed in unit");
  
  G4UIparameter* h3zValUnit = new G4UIparameter("zvalUnit", 's', true);
  h3zValUnit->SetGuidance("The unit applied to filled z-values and zvalMin0, zvalMax0");
  h3zValUnit->SetDefaultValue("none");
 
  G4UIparameter* h3zValFcn = new G4UIparameter("zvalFcn", 's', false);
  h3zValFcn->SetParameterCandidates("log log10 exp none");
  G4String fcnzGuidance = "The function applied to filled z-values (log, log10, exp, none).";
  h3zValFcn->SetGuidance(fcnzGuidance);
  h3zValFcn->SetDefaultValue("none");
 
  G4UIparameter* h3zValBinScheme = new G4UIparameter("zvalBinScheme", 's', true);
  G4String zbinSchemeGuidance = "The binning scheme (linear, log).";
  h3zValBinScheme->SetParameterCandidates("linear log");
  h3zValBinScheme->SetGuidance(zbinSchemeGuidance);
  h3zValBinScheme->SetDefaultValue("linear");
 
  fSetH3Cmd = new G4UIcommand("/analysis/h3/set", this);
  fSetH3Cmd->SetGuidance("Set parameters for the 3D histogram of #Id :");
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

//_____________________________________________________________________________
void G4H3Messenger::SetH3TitleCmd()
{
  G4UIparameter* h3Id = new G4UIparameter("idTitle", 'i', false);
  h3Id->SetGuidance("Histogram id");
  h3Id->SetParameterRange("idTitle>=0");

  G4UIparameter* h3Title = new G4UIparameter("h3Title", 's', true);
  h3Title->SetGuidance("Histogram title");
  h3Title->SetDefaultValue("none");

  fSetH3TitleCmd = new G4UIcommand("/analysis/h3/setTitle", this);
  fSetH3TitleCmd->SetGuidance("Set title for the 3D histogram of #Id");
  fSetH3TitleCmd->SetParameter(h3Id);
  fSetH3TitleCmd->SetParameter(h3Title);
  fSetH3TitleCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  

//_____________________________________________________________________________
void G4H3Messenger::SetH3XAxisCmd()
{
  G4UIparameter* h3Id = new G4UIparameter("idXaxis", 'i', false);
  h3Id->SetGuidance("Histogram id");
  h3Id->SetParameterRange("idXaxis>=0");

  G4UIparameter* h3XAxis = new G4UIparameter("h3Xaxis", 's', true);
  h3XAxis->SetGuidance("Histogram x-axis title");
  h3XAxis->SetDefaultValue("none");

  fSetH3XAxisCmd = new G4UIcommand("/analysis/h3/setXaxis", this);
  fSetH3XAxisCmd->SetGuidance("Set x-axis title for the 3D histogram of #Id");
  fSetH3XAxisCmd->SetParameter(h3Id);
  fSetH3XAxisCmd->SetParameter(h3XAxis);
  fSetH3XAxisCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  

//_____________________________________________________________________________
void G4H3Messenger::SetH3YAxisCmd()
{
  G4UIparameter* h3Id = new G4UIparameter("idYaxis", 'i', false);
  h3Id->SetGuidance("Histogram id");
  h3Id->SetParameterRange("idYaxis>=0");

  G4UIparameter* h3YAxis = new G4UIparameter("h3Yaxis", 's', true);
  h3YAxis->SetGuidance("Histogram y-axis title");
  h3YAxis->SetDefaultValue("none");

  fSetH3YAxisCmd = new G4UIcommand("/analysis/h3/setYaxis", this);
  fSetH3YAxisCmd->SetGuidance("Set y-axis title for the 3D histogram of #Id");
  fSetH3YAxisCmd->SetParameter(h3Id);
  fSetH3YAxisCmd->SetParameter(h3YAxis);
  fSetH3YAxisCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  

//_____________________________________________________________________________
void G4H3Messenger::SetH3ZAxisCmd()
{
  G4UIparameter* h3Id = new G4UIparameter("idYaxis", 'i', false);
  h3Id->SetGuidance("Histogram id");
  h3Id->SetParameterRange("idYaxis>=0");

  G4UIparameter* h3YAxis = new G4UIparameter("h3Yaxis", 's', true);
  h3YAxis->SetGuidance("Histogram y-axis title");
  h3YAxis->SetDefaultValue("none");

  fSetH3ZAxisCmd = new G4UIcommand("/analysis/h3/setZaxis", this);
  fSetH3ZAxisCmd->SetGuidance("Set y-axis title for the 3D histogram of #Id");
  fSetH3ZAxisCmd->SetParameter(h3Id);
  fSetH3ZAxisCmd->SetParameter(h3YAxis);
  fSetH3ZAxisCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
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
    Exception(command, parameters.size());
    return;
  }  

  if ( command == fCreateH3Cmd ) { 
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
    G4int znbins = G4UIcommand::ConvertToInt(parameters[counter++]); 
    G4double zvmin = G4UIcommand::ConvertToDouble(parameters[counter++]); 
    G4double zvmax = G4UIcommand::ConvertToDouble(parameters[counter++]); ; 
    G4String zsunit = parameters[counter++];
    G4String zsfcn = parameters[counter++];
    G4String zsbinScheme = parameters[counter++];
    G4double zunit = GetUnitValue(zsunit);
    fManager->CreateH3(name, title, 
                       xnbins, xvmin*xunit, xvmax*xunit,
                       ynbins, yvmin*yunit, yvmax*yunit, 
                       znbins, zvmin*zunit, zvmax*zunit, 
                       xsunit, ysunit, zsunit, xsfcn, ysfcn, zsfcn,
                       xsbinScheme, ysbinScheme, zsbinScheme);     
  }
  else if ( command == fSetH3Cmd ) {
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
    G4int znbins = G4UIcommand::ConvertToInt(parameters[counter++]); 
    G4double zvmin = G4UIcommand::ConvertToDouble(parameters[counter++]); 
    G4double zvmax = G4UIcommand::ConvertToDouble(parameters[counter++]); ; 
    G4String zsunit = parameters[counter++];
    G4String zsfcn = parameters[counter++];
    G4String zsbinScheme = parameters[counter++];
    G4double zunit = GetUnitValue(zsunit);
    fManager->SetH3(id, 
                    xnbins, xvmin*xunit, xvmax*xunit,
                    ynbins, yvmin*yunit, yvmax*yunit, 
                    znbins, zvmin*zunit, zvmax*zunit, 
                    xsunit, ysunit, zsunit, xsfcn, ysfcn, zsfcn,
                    xsbinScheme, ysbinScheme, zsbinScheme);     
  }
  else if ( command == fSetH3TitleCmd ) {
    G4int counter = 0;
    G4int id = G4UIcommand::ConvertToInt(parameters[counter++]); 
    G4String title = parameters[counter++];
    fManager->SetH3Title(id, title);     
  }
  else if ( command == fSetH3XAxisCmd ) {
    G4int counter = 0;
    G4int id = G4UIcommand::ConvertToInt(parameters[counter++]); 
    G4String xaxis = parameters[counter++];
    fManager->SetH3XAxisTitle(id, xaxis);     
  }
  else if ( command == fSetH3YAxisCmd ) {
    G4int counter = 0;
    G4int id = G4UIcommand::ConvertToInt(parameters[counter++]); 
    G4String yaxis = parameters[counter++];
    fManager->SetH3YAxisTitle(id, yaxis);     
  }
  else if ( command == fSetH3ZAxisCmd ) {
    G4int counter = 0;
    G4int id = G4UIcommand::ConvertToInt(parameters[counter++]); 
    G4String zaxis = parameters[counter++];
    fManager->SetH3ZAxisTitle(id, zaxis);     
  }
}  
