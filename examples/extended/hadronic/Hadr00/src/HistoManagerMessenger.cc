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
/// \file hadronic/Hadr00/src/HistoManagerMessenger.cc
/// \brief Implementation of the HistoManagerMessenger class
//
//
//
//
/////////////////////////////////////////////////////////////////////////
//
// HistoManagerMessenger
//
// Created: 20.06.08 V.Ivanchenko
//
// Modified:
//
////////////////////////////////////////////////////////////////////////
//

#include "HistoManagerMessenger.hh"

#include "HistoManager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "HistoManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManagerMessenger::HistoManagerMessenger(HistoManager * p)
:G4UImessenger(), fHisto(p)
{
  fbinCmd = new G4UIcmdWithAnInteger("/testhadr/nBinsE",this);
  fbinCmd->SetGuidance("Set number of bins for energy");
  fbinCmd->SetParameterName("NEbins",false);
  fbinCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fnOfAbsCmd = new G4UIcmdWithAnInteger("/testhadr/nBinsP",this);
  fnOfAbsCmd->SetGuidance("Set number of bins for momentum");
  fnOfAbsCmd->SetParameterName("NPbins",false);
  fnOfAbsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fpartCmd = new G4UIcmdWithAString("/testhadr/particle",this);
  fpartCmd->SetGuidance("Set particle name");
  fpartCmd->SetParameterName("Particle",false);
  fpartCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fcsCmd = new G4UIcmdWithAString("/testhadr/targetElm",this);
  fcsCmd->SetGuidance("Set element name");
  fcsCmd->SetParameterName("Elm",false);
  fcsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fe1Cmd = new G4UIcmdWithADoubleAndUnit("/testhadr/minEnergy",this);
  fe1Cmd->SetGuidance("Set min kinetic energy");
  fe1Cmd->SetParameterName("eMin",false);
  fe1Cmd->SetUnitCategory("Energy");
  fe1Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fe2Cmd = new G4UIcmdWithADoubleAndUnit("/testhadr/maxEnergy",this);
  fe2Cmd->SetGuidance("Set max kinetic energy");
  fe2Cmd->SetParameterName("eMax",false);
  fe2Cmd->SetUnitCategory("Energy");
  fe2Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fp1Cmd = new G4UIcmdWithADoubleAndUnit("/testhadr/minMomentum",this);
  fp1Cmd->SetGuidance("Set min momentum");
  fp1Cmd->SetParameterName("pMin",false);
  fp1Cmd->SetUnitCategory("Energy");
  fp1Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fp2Cmd = new G4UIcmdWithADoubleAndUnit("/testhadr/maxMomentum",this);
  fp2Cmd->SetGuidance("Set max momentum");
  fp2Cmd->SetParameterName("pMax",false);
  fp2Cmd->SetUnitCategory("Energy");
  fp2Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fverbCmd = new G4UIcmdWithAnInteger("/testhadr/verbose",this);
  fverbCmd->SetGuidance("Set verbose for ");
  fverbCmd->SetParameterName("verb",false);
  fverbCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fFCmd = new G4UIcmdWithAString("/testhadr/fileName",this);
  fFCmd->SetGuidance("set name for the histograms file");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManagerMessenger::~HistoManagerMessenger()
{
  delete fnOfAbsCmd;
  delete fpartCmd;
  delete fcsCmd;
  delete fe1Cmd;
  delete fe2Cmd;
  delete fp1Cmd;
  delete fp2Cmd;
  delete fverbCmd;
  delete fFCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManagerMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if( command == fbinCmd ) {
    fHisto->SetNumberOfBinsE(fbinCmd->GetNewIntValue(newValue));
  } else if( command == fnOfAbsCmd ) { 
    fHisto->SetNumberOfBinsP(fnOfAbsCmd->GetNewIntValue(newValue));
  } else if( command == fverbCmd ) {
    fHisto->SetVerbose(fverbCmd->GetNewIntValue(newValue));
  } else if( command == fpartCmd ) {
    fHisto->SetParticleName(newValue);
  } else if( command == fcsCmd ) {
    fHisto->SetElementName(newValue);
  } else if( command == fe1Cmd ) { 
    fHisto->SetMinKinEnergy(fe1Cmd->GetNewDoubleValue(newValue));
  } else if( command == fe2Cmd ) { 
    fHisto->SetMaxKinEnergy(fe2Cmd->GetNewDoubleValue(newValue));
  } else if( command == fp1Cmd ) { 
    fHisto->SetMinMomentum(fp1Cmd->GetNewDoubleValue(newValue));
  } else if( command == fp2Cmd ) { 
    fHisto->SetMaxMomentum(fp2Cmd->GetNewDoubleValue(newValue));
  } else if( command == fFCmd ) {
    fHisto->SetHistoName(newValue);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

