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
//
// $Id: DetectorMessenger.cc,v 1.1 2008-07-07 16:37:26 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
/////////////////////////////////////////////////////////////////////////
//
// DetectorMessenger
//
// Created: 20.06.08 V.Ivanchenko
//
// Modified:
//
////////////////////////////////////////////////////////////////////////
//

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "HistoManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:Detector(Det)
{
  testDir = new G4UIdirectory("/testhadr/");
  testDir->SetGuidance(" Hadronic Extended Example.");

  matCmd = new G4UIcmdWithAString("/testhadr/TargetMat",this);
  matCmd->SetGuidance("Select Material for the target");
  matCmd->SetParameterName("tMaterial",false);
  matCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  mat1Cmd = new G4UIcmdWithAString("/testhadr/WorldMat",this);
  mat1Cmd->SetGuidance("Select Material for world");
  mat1Cmd->SetParameterName("wMaterial",false);
  mat1Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  rCmd = new G4UIcmdWithADoubleAndUnit("/testhadr/TargetRadius",this);
  rCmd->SetGuidance("Set radius of the target");
  rCmd->SetParameterName("radius",false);
  rCmd->SetUnitCategory("Length");
  rCmd->SetRange("radius>0");
  rCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  lCmd = new G4UIcmdWithADoubleAndUnit("/testhadr/TargetLength",this);
  lCmd->SetGuidance("Set length of the target");
  lCmd->SetParameterName("length",false);
  lCmd->SetUnitCategory("Length");
  lCmd->SetRange("length>0");
  lCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  binCmd = new G4UIcmdWithAnInteger("/testhadr/nBinsE",this);
  binCmd->SetGuidance("Set number of bins for energy");
  binCmd->SetParameterName("NEbins",false);
  binCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  nOfAbsCmd = new G4UIcmdWithAnInteger("/testhadr/nBinsP",this);
  nOfAbsCmd->SetGuidance("Set number of bins for momentum");
  nOfAbsCmd->SetParameterName("NPbins",false);
  nOfAbsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  updateCmd = new G4UIcmdWithoutParameter("/testhadr/update",this);
  updateCmd->SetGuidance("Update geometry.");
  updateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  updateCmd->SetGuidance("if you changed geometrical value(s)");
  updateCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  partCmd = new G4UIcmdWithAString("/testhadr/particle",this);
  partCmd->SetGuidance("Set particle name");
  partCmd->SetParameterName("Particle",false);
  partCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  csCmd = new G4UIcmdWithAString("/testhadr/targetElm",this);
  csCmd->SetGuidance("Set element name");
  csCmd->SetParameterName("Elm",false);
  csCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  e1Cmd = new G4UIcmdWithADoubleAndUnit("/testhadr/minEnergy",this);
  e1Cmd->SetGuidance("Set min kinetic energy");
  e1Cmd->SetParameterName("eMin",false);
  e1Cmd->SetUnitCategory("Energy");
  e1Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  e2Cmd = new G4UIcmdWithADoubleAndUnit("/testhadr/maxEnergy",this);
  e2Cmd->SetGuidance("Set max kinetic energy");
  e2Cmd->SetParameterName("eMax",false);
  e2Cmd->SetUnitCategory("Energy");
  e2Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  p1Cmd = new G4UIcmdWithADoubleAndUnit("/testhadr/minMomentum",this);
  p1Cmd->SetGuidance("Set min momentum");
  p1Cmd->SetParameterName("pMin",false);
  p1Cmd->SetUnitCategory("Energy");
  p1Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  p2Cmd = new G4UIcmdWithADoubleAndUnit("/testhadr/maxMomentum",this);
  p2Cmd->SetGuidance("Set max momentum");
  p2Cmd->SetParameterName("pMax",false);
  p2Cmd->SetUnitCategory("Energy");
  p2Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  verbCmd = new G4UIcmdWithAnInteger("/testhadr/verbose",this);
  verbCmd->SetGuidance("Set verbose for ");
  verbCmd->SetParameterName("verb",false);
  verbCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete matCmd;
  delete mat1Cmd;
  delete rCmd;
  delete lCmd;
  delete nOfAbsCmd;
  delete updateCmd;
  delete testDir;
  delete partCmd;
  delete csCmd;
  delete e1Cmd;
  delete e2Cmd;
  delete p1Cmd;
  delete p2Cmd;
  delete verbCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  HistoManager* histo = HistoManager::GetPointer();
  if( command == matCmd ) {
    Detector->SetTargetMaterial(newValue);
  } else if( command == mat1Cmd ) {
    Detector->SetWorldMaterial(newValue);
  } else if( command == rCmd ) {
    Detector->SetTargetRadius(rCmd->GetNewDoubleValue(newValue));
  } else if( command == lCmd ) { 
    Detector->SetTargetLength(lCmd->GetNewDoubleValue(newValue));
  } else if( command == updateCmd ) {
    Detector->UpdateGeometry();
  } else if( command == binCmd ) {
    histo->SetNumberOfBinsE(binCmd->GetNewIntValue(newValue));
  } else if( command == nOfAbsCmd ) { 
    histo->SetNumberOfBinsP(nOfAbsCmd->GetNewIntValue(newValue));
  } else if( command == verbCmd ) {
    histo->SetVerbose(verbCmd->GetNewIntValue(newValue));
  } else if( command == partCmd ) {
    histo->SetParticleName(newValue);
  } else if( command == csCmd ) {
    histo->SetElementName(newValue);
  } else if( command == e1Cmd ) { 
    histo->SetMinKinEnergy(e1Cmd->GetNewDoubleValue(newValue));
  } else if( command == e2Cmd ) { 
    histo->SetMaxKinEnergy(e2Cmd->GetNewDoubleValue(newValue));
  } else if( command == p1Cmd ) { 
    histo->SetMinMomentum(p1Cmd->GetNewDoubleValue(newValue));
  } else if( command == p2Cmd ) { 
    histo->SetMaxMomentum(p2Cmd->GetNewDoubleValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

