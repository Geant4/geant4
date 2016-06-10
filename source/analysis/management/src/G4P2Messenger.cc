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

namespace {

void Exception(G4UIcommand* command, G4int nofParameters)
{
  G4ExceptionDescription description;
  description 
    << "Got wrong number of \"" << command->GetCommandName() 
    << "\" parameters: " << nofParameters
    << " instead of " << command->GetParameterEntries() 
    << " expected" << G4endl;
  G4Exception("G4P2Messenger::SetNewValue",
              "Analysis_W013", JustWarning, description);
}

}                  


//_____________________________________________________________________________
G4P2Messenger::G4P2Messenger(G4VAnalysisManager* manager)
  : G4UImessenger(),
    fManager(manager),
    fP2Dir(0),  
    fCreateP2Cmd(0),
    fSetP2Cmd(0),
    fSetP2TitleCmd(0), 
    fSetP2XAxisCmd(0), 
    fSetP2YAxisCmd(0)
{  
  fP2Dir = new G4UIdirectory("/analysis/p2/");
  fP2Dir->SetGuidance("2D profiles control");

  CreateP2Cmd();
  SetP2Cmd();
  
  SetP2TitleCmd();
  SetP2XAxisCmd();
  SetP2YAxisCmd();
  SetP2ZAxisCmd();
}

//_____________________________________________________________________________
G4P2Messenger::~G4P2Messenger()
{
  delete fCreateP2Cmd;
  delete fSetP2Cmd;
  delete fSetP2TitleCmd;  
  delete fSetP2XAxisCmd;  
  delete fSetP2YAxisCmd;  
  delete fSetP2ZAxisCmd;  
  delete fP2Dir;
}

//
// private functions
//

//_____________________________________________________________________________
void G4P2Messenger::CreateP2Cmd()
{
  G4UIparameter* p2Name = new G4UIparameter("name", 's', false);
  p2Name->SetGuidance("Profile name (label)");
  
  G4UIparameter* p2Title = new G4UIparameter("title", 's', false);
  p2Title->SetGuidance("Profile title");

  G4UIparameter* p2xNbins0 = new G4UIparameter("xnbins0", 'i', true);
  p2xNbins0->SetGuidance("Number of x-bins (default = 100)");
  p2xNbins0->SetGuidance("Can be reset with /analysis/p2/set command");
  p2xNbins0->SetDefaultValue(100);
  
  G4UIparameter* p2xValMin0 = new G4UIparameter("xvalMin0", 'd', true);
  p2xValMin0->SetGuidance("Minimum x-value, expressed in unit (default = 0.)");
  p2xValMin0->SetGuidance("Can be reset with /analysis/p2/set command");
  p2xValMin0->SetDefaultValue(0.);
  
  G4UIparameter* p2xValMax0 = new G4UIparameter("xvalMax0", 'd', true);
  p2xValMax0->SetGuidance("Maximum x-value, expressed in unit (default = 1.)");
  p2xValMax0->SetGuidance("Can be reset with /analysis/p2/set command");
  p2xValMax0->SetDefaultValue(1.);

  G4UIparameter* p2xValUnit0 = new G4UIparameter("xvalUnit0", 's', true);
  p2xValUnit0->SetGuidance("The unit applied to filled x-values and xvalMin0, xvalMax0");
  p2xValUnit0->SetDefaultValue("none");
  
  G4UIparameter* p2xValFcn0 = new G4UIparameter("xvalFcn0", 's', true);
  G4String fcnxGuidance = "The function applied to filled x-values (log, log10, exp, none).";
  p2xValFcn0->SetGuidance(fcnxGuidance);
  p2xValFcn0->SetParameterCandidates("log log10 exp none");
  p2xValFcn0->SetDefaultValue("none");
    
  G4UIparameter* p2xValBinScheme0 = new G4UIparameter("xvalBinScheme0", 's', true);
  G4String binSchemeGuidance = "The binning scheme (linear, log).";
  p2xValBinScheme0->SetParameterCandidates("linear log");
  p2xValBinScheme0->SetGuidance(binSchemeGuidance);
  p2xValBinScheme0->SetDefaultValue("linear");
  
  G4UIparameter* p2yNbins0 = new G4UIparameter("ynbins0", 'i', true);
  p2yNbins0->SetGuidance("Number of y-bins (default = 100)");
  p2yNbins0->SetGuidance("Can be reset with /analysis/p2/set command");
  p2yNbins0->SetDefaultValue(100);
  
  G4UIparameter* p2yValMin0 = new G4UIparameter("yvalMin0", 'd', true);
  p2yValMin0->SetGuidance("Minimum y-value, expressed in unit (default = 0.)");
  p2yValMin0->SetGuidance("Can be reset with /analysis/p2/set command");
  p2yValMin0->SetDefaultValue(0.);
  
  G4UIparameter* p2yValMax0 = new G4UIparameter("yvalMax0", 'd', true);
  p2yValMax0->SetGuidance("Maximum y-value, expressed in unit (default = 1.)");
  p2yValMax0->SetGuidance("Can be reset with /analysis/p2/set command");
  p2yValMax0->SetDefaultValue(1.);

  G4UIparameter* p2yValUnit0 = new G4UIparameter("yvalUnit0", 's', true);
  p2yValUnit0->SetGuidance("The unit applied to filled y-values and yvalMin0, yvalMax0");
  p2yValUnit0->SetDefaultValue("none");
  
  G4UIparameter* p2yValFcn0 = new G4UIparameter("yvalFcn0", 's', true);
  G4String fcnyGuidance = "The function applied to filled y-values (log, log10, exp, none).";
  p2yValFcn0->SetGuidance(fcnyGuidance);
  p2yValFcn0->SetParameterCandidates("log log10 exp none");
  p2yValFcn0->SetDefaultValue("none");
    
  G4UIparameter* p2yValBinScheme0 = new G4UIparameter("yvalBinScheme0", 's', true);
  p2yValBinScheme0->SetParameterCandidates("linear log");
  p2yValBinScheme0->SetGuidance(binSchemeGuidance);
  p2yValBinScheme0->SetDefaultValue("linear");
  
  G4UIparameter* p2zValMin0 = new G4UIparameter("zvalMin0", 'd', true);
  p2zValMin0->SetGuidance("Minimum z-value, expressed in unit (default = 0.)");
  p2zValMin0->SetGuidance("Can be reset with /analysis/p2/set command");
  p2zValMin0->SetDefaultValue(0.);
  
  G4UIparameter* p2zValMax0 = new G4UIparameter("zvalMax0", 'd', true);
  p2zValMax0->SetGuidance("Maximum z-value, expressed in unit (default = 1.)");
  p2zValMax0->SetGuidance("Can be reset with /analysis/p2/set command");
  p2zValMax0->SetDefaultValue(1.);

  G4UIparameter* p2zValUnit0 = new G4UIparameter("zvalUnit0", 's', true);
  p2zValUnit0->SetGuidance("The unit applied to filled z-values and zvalMin0, zvalMax0");
  p2zValUnit0->SetDefaultValue("none");
  
  G4UIparameter* p2zValFcn0 = new G4UIparameter("zvalFcn0", 's', true);
  G4String fcnzGuidance = "The function applied to filled z-values (log, log10, exp, none).";
  p2zValFcn0->SetGuidance(fcnzGuidance);
  p2zValFcn0->SetParameterCandidates("log log10 exp none");
  p2zValFcn0->SetDefaultValue("none");
  
  fCreateP2Cmd = new G4UIcommand("/analysis/p2/create", this);
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
  G4UIparameter* p2Id = new G4UIparameter("id", 'i', false);
  p2Id->SetGuidance("Profile id");
  p2Id->SetParameterRange("id>=0");
  
  G4UIparameter* p2xNbins = new G4UIparameter("xnbins", 'i', false);
  p2xNbins->SetGuidance("Number of x-bins");
  
  G4UIparameter* p2xValMin = new G4UIparameter("xvalMin", 'd', false);
  p2xValMin->SetGuidance("Minimum x-value, expressed in unit");
  
  G4UIparameter* p2xValMax = new G4UIparameter("xvalMax", 'd', false);
  p2xValMax->SetGuidance("Maximum x-value, expressed in unit");
  
  G4UIparameter* p2xValUnit = new G4UIparameter("xvalUnit", 's', false);
  p2xValUnit->SetGuidance("The unit applied to filled x-values and xvalMin0, xvalMax0");
  p2xValUnit->SetDefaultValue("none");
 
  G4UIparameter* p2xValFcn = new G4UIparameter("xvalFcn", 's', false);
  p2xValFcn->SetParameterCandidates("log log10 exp none");
  G4String fcnxGuidance = "The function applied to filled x-values (log, log10, exp, none).";
  p2xValFcn->SetGuidance(fcnxGuidance);
  p2xValFcn->SetDefaultValue("none");
    
  G4UIparameter* p2xValBinScheme = new G4UIparameter("xvalBinScheme", 's', true);
  G4String binSchemeGuidance = "The binning scheme (linear, log).";
  p2xValBinScheme->SetParameterCandidates("linear log");
  p2xValBinScheme->SetGuidance(binSchemeGuidance);
  p2xValBinScheme->SetDefaultValue("linear");
 
  G4UIparameter* p2yNbins = new G4UIparameter("nybins", 'i', false);
  p2yNbins->SetGuidance("Number of y-bins");
  
  G4UIparameter* p2yValMin = new G4UIparameter("yvalMin", 'd', false);
  p2yValMin->SetGuidance("Minimum y-value, expressed in unit");
  
  G4UIparameter* p2yValMax = new G4UIparameter("yvalMax", 'd', false);
  p2yValMax->SetGuidance("Maximum y-value, expressed in unit");
  
  G4UIparameter* p2yValUnit = new G4UIparameter("yvalUnit", 's', true);
  p2yValUnit->SetGuidance("The unit applied to filled y-values and yvalMin0, yvalMax0");
  p2yValUnit->SetDefaultValue("none");
 
  G4UIparameter* p2yValFcn = new G4UIparameter("yvalFcn", 's', false);
  p2yValFcn->SetParameterCandidates("log log10 exp none");
  G4String fcnyGuidance = "The function applied to filled y-values (log, log10, exp, none).";
  p2yValFcn->SetGuidance(fcnyGuidance);
  p2yValFcn->SetDefaultValue("none");
    
  G4UIparameter* p2yValBinScheme = new G4UIparameter("yvalBinScheme", 's', true);
  p2yValBinScheme->SetParameterCandidates("linear log");
  p2yValBinScheme->SetGuidance(binSchemeGuidance);
  p2yValBinScheme->SetDefaultValue("linear");
 
  G4UIparameter* p2zValMin = new G4UIparameter("zvalMin", 'd', false);
  p2zValMin->SetGuidance("Minimum z-value, expressed in unit");
  
  G4UIparameter* p2zValMax = new G4UIparameter("zvalMax", 'd', false);
  p2zValMax->SetGuidance("Maximum z-value, expressed in unit");
  
  G4UIparameter* p2zValUnit = new G4UIparameter("zvalUnit", 's', true);
  p2zValUnit->SetGuidance("The unit applied to filled z-values and zvalMin0, zvalMax0");
  p2zValUnit->SetDefaultValue("none");
 
  G4UIparameter* p2zValFcn = new G4UIparameter("zvalFcn", 's', false);
  p2zValFcn->SetParameterCandidates("log log10 exp none");
  G4String fcnzGuidance = "The function applied to filled z-values (log, log10, exp, none).";
  p2zValFcn->SetGuidance(fcnzGuidance);
  p2zValFcn->SetDefaultValue("none");
 
  fSetP2Cmd = new G4UIcommand("/analysis/p2/set", this);
  fSetP2Cmd->SetGuidance("Set parameters for the 2D profile of #Id :");
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

//_____________________________________________________________________________
void G4P2Messenger::SetP2TitleCmd()
{
  G4UIparameter* p2Id = new G4UIparameter("idTitle", 'i', false);
  p2Id->SetGuidance("Profile id");
  p2Id->SetParameterRange("idTitle>=0");

  G4UIparameter* p2Title = new G4UIparameter("p2Title", 's', true);
  p2Title->SetGuidance("Profile title");
  p2Title->SetDefaultValue("none");

  fSetP2TitleCmd = new G4UIcommand("/analysis/p2/setTitle", this);
  fSetP2TitleCmd->SetGuidance("Set title for the 2D profile of #Id");
  fSetP2TitleCmd->SetParameter(p2Id);
  fSetP2TitleCmd->SetParameter(p2Title);
  fSetP2TitleCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  

//_____________________________________________________________________________
void G4P2Messenger::SetP2XAxisCmd()
{
  G4UIparameter* p2Id = new G4UIparameter("idXaxis", 'i', false);
  p2Id->SetGuidance("Profile id");
  p2Id->SetParameterRange("idXaxis>=0");

  G4UIparameter* p2XAxis = new G4UIparameter("p2Xaxis", 's', true);
  p2XAxis->SetGuidance("Profile x-axis title");
  p2XAxis->SetDefaultValue("none");

  fSetP2XAxisCmd = new G4UIcommand("/analysis/p2/setXaxis", this);
  fSetP2XAxisCmd->SetGuidance("Set x-axis title for the 2D profile of #Id");
  fSetP2XAxisCmd->SetParameter(p2Id);
  fSetP2XAxisCmd->SetParameter(p2XAxis);
  fSetP2XAxisCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  

//_____________________________________________________________________________
void G4P2Messenger::SetP2YAxisCmd()
{
  G4UIparameter* p2Id = new G4UIparameter("idYaxis", 'i', false);
  p2Id->SetGuidance("Profile id");
  p2Id->SetParameterRange("idYaxis>=0");

  G4UIparameter* p2YAxis = new G4UIparameter("p2Yaxis", 's', true);
  p2YAxis->SetGuidance("Profile y-axis title");
  p2YAxis->SetDefaultValue("none");

  fSetP2YAxisCmd = new G4UIcommand("/analysis/p2/setYaxis", this);
  fSetP2YAxisCmd->SetGuidance("Set y-axis title for the 2D profile of #Id");
  fSetP2YAxisCmd->SetParameter(p2Id);
  fSetP2YAxisCmd->SetParameter(p2YAxis);
  fSetP2YAxisCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  

//_____________________________________________________________________________
void G4P2Messenger::SetP2ZAxisCmd()
{
  G4UIparameter* p2Id = new G4UIparameter("idZaxis", 'i', false);
  p2Id->SetGuidance("Profile id");
  p2Id->SetParameterRange("idZaxis>=0");

  G4UIparameter* p2ZAxis = new G4UIparameter("p2Zaxis", 's', true);
  p2ZAxis->SetGuidance("Profile z-axis title");
  p2ZAxis->SetDefaultValue("none");

  fSetP2ZAxisCmd = new G4UIcommand("/analysis/p2/setZaxis", this);
  fSetP2ZAxisCmd->SetParameter(p2Id);
  fSetP2ZAxisCmd->SetParameter(p2ZAxis);
  fSetP2ZAxisCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
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
    Exception(command, parameters.size());
    return;
  }  

  if ( command == fCreateP2Cmd ) { 
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
    G4double zvmin = G4UIcommand::ConvertToDouble(parameters[counter++]); 
    G4double zvmax = G4UIcommand::ConvertToDouble(parameters[counter++]); ; 
    G4String zsunit = parameters[counter++];
    G4String zsfcn = parameters[counter++];
    G4double zunit = GetUnitValue(zsunit);
    fManager->CreateP2(name, title, 
                       xnbins, xvmin*xunit, xvmax*xunit,
                       ynbins, yvmin*yunit, yvmax*yunit, 
                       zvmin*zunit, zvmax*zunit, 
                       xsunit, ysunit, zsunit, xsfcn, ysfcn, zsfcn,
                       xsbinScheme, ysbinScheme);     
  }
  else if ( command == fSetP2Cmd ) {
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
    G4double zvmin = G4UIcommand::ConvertToDouble(parameters[counter++]); 
    G4double zvmax = G4UIcommand::ConvertToDouble(parameters[counter++]); ; 
    G4String zsunit = parameters[counter++];
    G4String zsfcn = parameters[counter++];
    G4double zunit = GetUnitValue(zsunit);
    fManager->SetP2(id, 
                    xnbins, xvmin*xunit, xvmax*xunit,
                    ynbins, yvmin*yunit, yvmax*yunit, 
                    zvmin*zunit, zvmax*zunit, 
                    xsunit, ysunit, zsunit, xsfcn, ysfcn, zsfcn,     
                    xsbinScheme, ysbinScheme);     
  }
  else if ( command == fSetP2TitleCmd ) {
    G4int counter = 0;
    G4int id = G4UIcommand::ConvertToInt(parameters[counter++]); 
    G4String title = parameters[counter++];
    fManager->SetP2Title(id, title);     
  }
  else if ( command == fSetP2XAxisCmd ) {
    G4int counter = 0;
    G4int id = G4UIcommand::ConvertToInt(parameters[counter++]); 
    G4String xaxis = parameters[counter++];
    fManager->SetP2XAxisTitle(id, xaxis);     
  }
  else if ( command == fSetP2YAxisCmd ) {
    G4int counter = 0;
    G4int id = G4UIcommand::ConvertToInt(parameters[counter++]); 
    G4String yaxis = parameters[counter++];
    fManager->SetP2YAxisTitle(id, yaxis);     
  }
  else if ( command == fSetP2ZAxisCmd ) {
    G4int counter = 0;
    G4int id = G4UIcommand::ConvertToInt(parameters[counter++]); 
    G4String zaxis = parameters[counter++];
    fManager->SetP2ZAxisTitle(id, zaxis);     
  }
}  
