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

// Author: Ivana Hrivnacova, 04/07/2012  (ivana@ipno.in2p3.fr)
//
// This messenger class is a generalization of the HistoMessenger class,
// originally developed for the extended/electromagnetic examples
// by Michel Maire (michel.maire@lapp.in2p3.fr)

#include "G4AnalysisMessenger.hh"
#include "G4VAnalysisManager.hh"
#include "G4HnInformation.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UnitsTable.hh"

#include <iostream>

//_____________________________________________________________________________
G4AnalysisMessenger::G4AnalysisMessenger(G4VAnalysisManager* manager)
  : G4UImessenger(),
    fManager(manager),
    fAnalysisDir(0),  
    fSetFileNameCmd(0),
    fSetHistoDirNameCmd(0),
    fSetNtupleDirNameCmd(0),
    fVerboseCmd(0),
    fH1Dir(0),  
    fCreateH1Cmd(0),
    fSetH1Cmd(0),
    fSetH1Ascii(0) 
{  
  fAnalysisDir = new G4UIdirectory("/analysis/");
  fAnalysisDir->SetGuidance("analysis control");

  fSetFileNameCmd = new G4UIcmdWithAString("/analysis/setFileName",this);
  fSetFileNameCmd->SetGuidance("Set name for the histograms & ntuple file");
  fSetFileNameCmd->SetParameterName("Filename", false);
  fSetFileNameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  fSetHistoDirNameCmd = new G4UIcmdWithAString("/analysis/setHistoDirName",this);
  fSetHistoDirNameCmd->SetGuidance("Set name for the histograms directory");
  fSetHistoDirNameCmd->SetParameterName("HistoDirName", false);
  fSetHistoDirNameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  fSetNtupleDirNameCmd = new G4UIcmdWithAString("/analysis/setNtupleDirName",this);
  fSetNtupleDirNameCmd->SetGuidance("Set name for the ntuple directory");
  fSetNtupleDirNameCmd->SetParameterName("NtupleDirName", false);
  fSetNtupleDirNameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  fVerboseCmd = new G4UIcmdWithAnInteger("/analysis/verbose",this);
  fVerboseCmd->SetGuidance("Set verbose level");
  fVerboseCmd->SetParameterName("VerboseLevel",false);
  fVerboseCmd->SetRange("VerboseLevel>=0 && VerboseLevel<=3");

  fH1Dir = new G4UIdirectory("/analysis/h1/");
  fH1Dir->SetGuidance("1D histograms control");

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
  h1ValUnit0->SetGuidance("The unit of valMin0 and valMax0");
  h1ValUnit0->SetDefaultValue("none");
  
  fCreateH1Cmd = new G4UIcommand("/analysis/h1/create", this);
  fCreateH1Cmd->SetGuidance("Create 1D histogram");
  fCreateH1Cmd->SetParameter(h1Name);
  fCreateH1Cmd->SetParameter(h1Title);
  fCreateH1Cmd->SetParameter(h1Nbins0);
  fCreateH1Cmd->SetParameter(h1ValMin0);
  fCreateH1Cmd->SetParameter(h1ValMax0);
  fCreateH1Cmd->SetParameter(h1ValUnit0);
  fCreateH1Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
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
  h1ValUnit->SetGuidance("The unit of valMin and valMax");
  h1ValUnit->SetDefaultValue("none");
 
  fSetH1Cmd = new G4UIcommand("/analysis/h1/set", this);
  fSetH1Cmd->SetGuidance("Set parameters for the 1D histogram of #Id :");
  fSetH1Cmd->SetGuidance("  nbins; valMin; valMax; unit (of vmin and vmax)");
  fSetH1Cmd->SetParameter(h1Id);
  fSetH1Cmd->SetParameter(h1Nbins);
  fSetH1Cmd->SetParameter(h1ValMin);
  fSetH1Cmd->SetParameter(h1ValMax);
  fSetH1Cmd->SetParameter(h1ValUnit);
  fSetH1Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fSetH1Ascii = new G4UIcmdWithAnInteger("/analysis/h1/setAscii",this);
  fSetH1Ascii->SetGuidance("Print 1D histogram of #Id on ascii file.");
  fSetH1Ascii->SetParameterName("Id",false);
  fSetH1Ascii->SetRange("Id>=0");
  fSetH1Ascii->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//_____________________________________________________________________________
G4AnalysisMessenger::~G4AnalysisMessenger()
{
  delete fSetFileNameCmd;
  delete fSetHistoDirNameCmd;
  delete fSetNtupleDirNameCmd;
  delete fVerboseCmd;
  delete fCreateH1Cmd;
  delete fSetH1Cmd;
  delete fSetH1Ascii;  
  delete fH1Dir;
  delete fAnalysisDir;
}

//_____________________________________________________________________________
void G4AnalysisMessenger::SetNewValue(G4UIcommand* command, G4String newValues)
{
  if ( command == fSetFileNameCmd ) {
    G4cout << "Set file name: " << newValues << G4endl;
    fManager->SetFileName(newValues);
  }  
  else if ( command == fSetHistoDirNameCmd ) {
    fManager->SetHistoDirectoryName(newValues);
  }  
  else if ( command == fSetNtupleDirNameCmd ) {
    fManager->SetNtupleDirectoryName(newValues);
  }  
  else if ( command == fVerboseCmd ) {
    fManager->SetVerboseLevel(fVerboseCmd->GetNewIntValue(newValues));
  }  
  else if ( command == fCreateH1Cmd ) { 
    G4String name, title;
    G4int nbins; 
    G4double vmin,vmax; 
    G4String sunit;
    std::istringstream is(newValues.data());
    is >> name >> title >> nbins >> vmin >> vmax >> sunit;

    G4double unit = 1.;
    if ( sunit != "none" ) {
      unit = G4UnitDefinition::GetValueOf(sunit);
    }  
    fManager->CreateH1(name, title, nbins, vmin, vmax, unit);     
  }
  else if ( command == fSetH1Cmd ) {
    G4int id; 
    G4int nbins; 
    G4double vmin, vmax; 
    G4String sunit;
    std::istringstream is(newValues.data());
    is >> id >> nbins >> vmin >> vmax >> sunit;
    G4double unit = 1.;
    if ( sunit != "none" ) unit = G4UnitDefinition::GetValueOf(sunit);
    if ( unit == 0. ) unit = 1.;  
    fManager->SetH1(id, nbins, vmin, vmax, unit);     
  }
  else if ( command == fSetH1Ascii ) {
    G4int id = fSetH1Ascii->GetNewIntValue(newValues);
    fManager->SetAscii(G4VAnalysisManager::kH1, id, true); 
  }      
}  
