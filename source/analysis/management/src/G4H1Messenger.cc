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
// $Id: G4H1Messenger.cc 66310 2012-12-17 11:56:35Z ihrivnac $

// Author: Ivana Hrivnacova, 24/06/2013  (ivana@ipno.in2p3.fr)
//
// This messenger class is a generalization of the HistoMessenger class,
// originally developed for the extended/electromagnetic examples
// by Michel Maire (michel.maire@lapp.in2p3.fr)

#include "G4H1Messenger.hh"
#include "G4VAnalysisManager.hh"
#include "G4AnalysisUtilities.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4Tokenizer.hh"

#include <iostream>
#include <vector>

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
  G4Exception("G4H1Messenger::SetNewValue",
              "Analysis_W013", JustWarning, description);
}

}                  

//_____________________________________________________________________________
G4H1Messenger::G4H1Messenger(G4VAnalysisManager* manager)
  : G4UImessenger(),
    fManager(manager),
    fH1Dir(0),  
    fCreateH1Cmd(0),
    fSetH1Cmd(0),
    fSetH1TitleCmd(0), 
    fSetH1XAxisCmd(0), 
    fSetH1YAxisCmd(0)
{  
  fH1Dir = new G4UIdirectory("/analysis/h1/");
  fH1Dir->SetGuidance("1D histograms control");

  CreateH1Cmd();
  SetH1Cmd();
  SetH1TitleCmd();
  SetH1XAxisCmd();
  SetH1YAxisCmd();
}

//_____________________________________________________________________________
G4H1Messenger::~G4H1Messenger()
{
  delete fCreateH1Cmd;
  delete fSetH1Cmd;
  delete fSetH1TitleCmd;  
  delete fSetH1XAxisCmd;  
  delete fSetH1YAxisCmd;  
  delete fH1Dir;
}

//
// private functions
//

//_____________________________________________________________________________
void G4H1Messenger::CreateH1Cmd()
{
  G4UIparameter* h1Name = new G4UIparameter("name", 's', false);
  h1Name->SetGuidance("Histogram name (label)");
  
  G4UIparameter* h1Title = new G4UIparameter("title", 's', false);
  h1Title->SetGuidance("Histogram title");

  G4UIparameter* h1Nbins0 = new G4UIparameter("nbins0", 'i', true);
  h1Nbins0->SetGuidance("Number of bins (default = 100)");
  h1Nbins0->SetGuidance("Can be reset with /analysis/h1/set command");
  h1Nbins0->SetDefaultValue(100);
  
  G4UIparameter* h1ValMin0 = new G4UIparameter("valMin0", 'd', true);
  h1ValMin0->SetGuidance("Minimum value, expressed in unit (default = 0.)");
  h1ValMin0->SetGuidance("Can be reset with /analysis/h1/set command");
  h1ValMin0->SetDefaultValue(0.);
  
  G4UIparameter* h1ValMax0 = new G4UIparameter("valMax0", 'd', true);
  h1ValMax0->SetGuidance("Maximum value, expressed in unit (default = 1.)");
  h1ValMax0->SetGuidance("Can be reset with /analysis/h1/set command");
  h1ValMax0->SetDefaultValue(1.);

  G4UIparameter* h1ValUnit0 = new G4UIparameter("valUnit0", 's', true);
  h1ValUnit0->SetGuidance("The unit applied to filled values and valMin0, valMax0");
  h1ValUnit0->SetDefaultValue("none");
  
  G4UIparameter* h1ValFcn0 = new G4UIparameter("valFcn0", 's', true);
  G4String fcnGuidance = "The function applied to filled values (log, log10, exp).\n";
  fcnGuidance += "Note that the unit parameter cannot be omitted in this case,\n";
  fcnGuidance += "but none value should be used instead.";
  h1ValFcn0->SetGuidance(fcnGuidance);
  h1ValFcn0->SetParameterCandidates("log log10 exp none");
  h1ValFcn0->SetDefaultValue("none");
  
  G4UIparameter* h1ValBinScheme0 = new G4UIparameter("valBinScheme0", 's', true);
  G4String binSchemeGuidance = "The binning scheme (linear, log).\n";
  h1ValBinScheme0->SetParameterCandidates("linear log");
  binSchemeGuidance 
    += "Note that the unit and fcn parameters cannot be omitted in this case,\n";
  binSchemeGuidance += "but none value should be used instead.";
  h1ValBinScheme0->SetGuidance(binSchemeGuidance);
  h1ValBinScheme0->SetDefaultValue("linear");
  
  fCreateH1Cmd = new G4UIcommand("/analysis/h1/create", this);
  fCreateH1Cmd->SetGuidance("Create 1D histogram");
  fCreateH1Cmd->SetParameter(h1Name);
  fCreateH1Cmd->SetParameter(h1Title);
  fCreateH1Cmd->SetParameter(h1Nbins0);
  fCreateH1Cmd->SetParameter(h1ValMin0);
  fCreateH1Cmd->SetParameter(h1ValMax0);
  fCreateH1Cmd->SetParameter(h1ValUnit0);
  fCreateH1Cmd->SetParameter(h1ValFcn0);
  fCreateH1Cmd->SetParameter(h1ValBinScheme0);
  fCreateH1Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  


//_____________________________________________________________________________
void G4H1Messenger::SetH1Cmd()
{
  G4UIparameter* h1Id = new G4UIparameter("id", 'i', false);
  h1Id->SetGuidance("Histogram id");
  h1Id->SetParameterRange("id>=0");
  
  G4UIparameter* h1Nbins = new G4UIparameter("nbins", 'i', false);
  h1Nbins->SetGuidance("Number of bins");
  
  G4UIparameter* h1ValMin = new G4UIparameter("valMin", 'd', false);
  h1ValMin->SetGuidance("Minimum value, expressed in unit");
  
  G4UIparameter* h1ValMax = new G4UIparameter("valMax", 'd', false);
  h1ValMax->SetGuidance("Maximum value, expressed in unit");
  
  G4UIparameter* h1ValUnit = new G4UIparameter("valUnit", 's', true);
  h1ValUnit->SetGuidance("The unit applied to filled values and valMin, valMax");
  h1ValUnit->SetDefaultValue("none");

  G4UIparameter* h1ValFcn = new G4UIparameter("valFcn", 's', true);
  h1ValFcn->SetParameterCandidates("log log10 exp none");
  G4String fcnGuidance = "The function applied to filled values (log, log10, exp, none).\n";
  fcnGuidance += "Note that the unit parameter cannot be omitted in this case,\n";
  fcnGuidance += "but none value should be used instead.";
  h1ValFcn->SetGuidance(fcnGuidance);
  h1ValFcn->SetDefaultValue("none");

  G4UIparameter* h1ValBinScheme = new G4UIparameter("valBinScheme", 's', true);
  h1ValBinScheme->SetParameterCandidates("linear log");
  G4String binSchemeGuidance = "The binning scheme (linear, log).\n";
  binSchemeGuidance 
    += "Note that the unit and fcn parameters cannot be omitted in this case,\n";
  binSchemeGuidance += "but none value should be used instead.";
  h1ValBinScheme->SetGuidance(binSchemeGuidance);
  h1ValBinScheme->SetDefaultValue("linear");
 
  fSetH1Cmd = new G4UIcommand("/analysis/h1/set", this);
  fSetH1Cmd->SetGuidance("Set parameters for the 1D histogram of #Id :");
  fSetH1Cmd->SetGuidance("  nbins; valMin; valMax; unit; function; binScheme");
  fSetH1Cmd->SetParameter(h1Id);
  fSetH1Cmd->SetParameter(h1Nbins);
  fSetH1Cmd->SetParameter(h1ValMin);
  fSetH1Cmd->SetParameter(h1ValMax);
  fSetH1Cmd->SetParameter(h1ValUnit);
  fSetH1Cmd->SetParameter(h1ValFcn);
  fSetH1Cmd->SetParameter(h1ValBinScheme);
  fSetH1Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  

//_____________________________________________________________________________
void G4H1Messenger::SetH1TitleCmd()
{
  G4UIparameter* h1Id = new G4UIparameter("idTitle", 'i', false);
  h1Id->SetGuidance("Histogram id");
  h1Id->SetParameterRange("idTitle>=0");

  G4UIparameter* h1Title = new G4UIparameter("h1Title", 's', true);
  h1Title->SetGuidance("Histogram title");
  h1Title->SetDefaultValue("none");

  fSetH1TitleCmd = new G4UIcommand("/analysis/h1/setTitle", this);
  fSetH1TitleCmd->SetGuidance("Set title for the 1D histogram of #Id");
  fSetH1TitleCmd->SetParameter(h1Id);
  fSetH1TitleCmd->SetParameter(h1Title);
  fSetH1TitleCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  

//_____________________________________________________________________________
void G4H1Messenger::SetH1XAxisCmd()
{
  G4UIparameter* h1Id = new G4UIparameter("idXaxis", 'i', false);
  h1Id->SetGuidance("Histogram id");
  h1Id->SetParameterRange("idXaxis>=0");

  G4UIparameter* h1XAxis = new G4UIparameter("h1Xaxis", 's', true);
  h1XAxis->SetGuidance("Histogram x-axis title");
  h1XAxis->SetDefaultValue("none");

  fSetH1XAxisCmd = new G4UIcommand("/analysis/h1/setXaxis", this);
  fSetH1XAxisCmd->SetGuidance("Set x-axis title for the 1D histogram of #Id");
  fSetH1XAxisCmd->SetParameter(h1Id);
  fSetH1XAxisCmd->SetParameter(h1XAxis);
  fSetH1XAxisCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  

//_____________________________________________________________________________
void G4H1Messenger::SetH1YAxisCmd()
{
  G4UIparameter* h1Id = new G4UIparameter("idYaxis", 'i', false);
  h1Id->SetGuidance("Histogram id");
  h1Id->SetParameterRange("idYaxis>=0");

  G4UIparameter* h1YAxis = new G4UIparameter("h1Yaxis", 's', true);
  h1YAxis->SetGuidance("Histogram y-axis title");
  h1YAxis->SetDefaultValue("none");

  fSetH1YAxisCmd = new G4UIcommand("/analysis/h1/setYaxis", this);
  fSetH1YAxisCmd->SetGuidance("Set y-axis title for the 1D histogram of #Id");
  fSetH1YAxisCmd->SetParameter(h1Id);
  fSetH1YAxisCmd->SetParameter(h1YAxis);
  fSetH1YAxisCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}  

//
// public functions
//

//_____________________________________________________________________________
void G4H1Messenger::SetNewValue(G4UIcommand* command, G4String newValues)
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

  if ( command == fCreateH1Cmd ) { 
    G4int counter = 0;
    G4String name = parameters[counter++];
    G4String title = parameters[counter++];
    G4int nbins = G4UIcommand::ConvertToInt(parameters[counter++]); 
    G4double vmin = G4UIcommand::ConvertToDouble(parameters[counter++]); 
    G4double vmax = G4UIcommand::ConvertToDouble(parameters[counter++]); ; 
    G4String sunit = parameters[counter++];
    G4String sfcn = parameters[counter++];
    G4String sbinScheme = parameters[counter++];
    G4double unit = GetUnitValue(sunit);
    fManager->CreateH1(name, title, nbins, vmin*unit, vmax*unit, sunit, sfcn, sbinScheme); 
  }
  else if ( command == fSetH1Cmd ) {
    G4int counter = 0;
    G4int id = G4UIcommand::ConvertToInt(parameters[counter++]);
    G4int nbins = G4UIcommand::ConvertToInt(parameters[counter++]); 
    G4double vmin = G4UIcommand::ConvertToDouble(parameters[counter++]); 
    G4double vmax = G4UIcommand::ConvertToDouble(parameters[counter++]); 
    G4String sunit = parameters[counter++];
    G4String sfcn = parameters[counter++];
    G4String sbinScheme = parameters[counter++];
    G4double unit = GetUnitValue(sunit);
    fManager->SetH1(id, nbins, vmin*unit, vmax*unit, sunit, sfcn, sbinScheme); 
  }
  else if ( command == fSetH1TitleCmd ) {
    G4int counter = 0;
    G4int id = G4UIcommand::ConvertToInt(parameters[counter++]); 
    G4String title = parameters[counter++];
    fManager->SetH1Title(id, title);     
  }
  else if ( command == fSetH1XAxisCmd ) {
    G4int counter = 0;
    G4int id = G4UIcommand::ConvertToInt(parameters[counter++]); 
    G4String xaxis = parameters[counter++];
    fManager->SetH1XAxisTitle(id, xaxis);     
  }
  else if ( command == fSetH1YAxisCmd ) {
    G4int counter = 0;
    G4int id = G4UIcommand::ConvertToInt(parameters[counter++]); 
    G4String yaxis = parameters[counter++];
    fManager->SetH1YAxisTitle(id, yaxis);     
  }
}  
