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
/// \file hadronic/Hadr02/src/DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
// $Id$
//
/////////////////////////////////////////////////////////////////////////
//
// DetectorMessenger
//
// Created: 31.01.03 V.Ivanchenko
//
// Modified:
// 04.06.2006 Adoptation of hadr01 (V.Ivanchenko)
// 16.11.2006 Add beamCmd (V.Ivanchenko)
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

  ionCmd = new G4UIcmdWithAString("/testhadr/ionPhysics",this);
  ionCmd->SetGuidance("Select ion Physics");
  ionCmd->SetParameterName("DPMJET",false);
  ionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

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

  binCmd = new G4UIcmdWithAnInteger("/testhadr/NumberOfBinsE",this);
  binCmd->SetGuidance("Set number of bins for Energy");
  binCmd->SetParameterName("NEbins",false);
  binCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  nOfAbsCmd = new G4UIcmdWithAnInteger("/testhadr/NumberDivZ",this);
  nOfAbsCmd->SetGuidance("Set number of slices");
  nOfAbsCmd->SetParameterName("NZ",false);
  nOfAbsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  edepCmd = new G4UIcmdWithADoubleAndUnit("/testhadr/MaxEdep",this);
  edepCmd->SetGuidance("Set max energy in histogram");
  edepCmd->SetParameterName("edep",false);
  edepCmd->SetUnitCategory("Energy");
  edepCmd->SetRange("edep>0");
  edepCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  beamCmd = new G4UIcmdWithABool("/testhadr/DefaultBeamPosition",this);
  beamCmd->SetGuidance("show inelastic and elastic cross sections");

  verbCmd = new G4UIcmdWithAnInteger("/testhadr/Verbose",this);
  verbCmd->SetGuidance("Set verbose for ");
  verbCmd->SetParameterName("verb",false);
  verbCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete matCmd;
  delete mat1Cmd;
  delete ionCmd;
  delete rCmd;
  delete lCmd;
  delete nOfAbsCmd;
  delete testDir;
  delete beamCmd;
  delete verbCmd;
  delete edepCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  HistoManager* h = HistoManager::GetPointer();
  if( command == matCmd ) {
    Detector->SetTargetMaterial(newValue);
  } else if( command == mat1Cmd ) {
    Detector->SetWorldMaterial(newValue);
  } else if( command == ionCmd ) {
    h->SetIonPhysics(newValue);
  } else if( command == rCmd ) {
    Detector->SetTargetRadius(rCmd->GetNewDoubleValue(newValue));
  } else if( command == lCmd ) {
    h->SetTargetLength(lCmd->GetNewDoubleValue(newValue));
  } else if( command == nOfAbsCmd ) {
    h->SetNumberOfSlices(nOfAbsCmd->GetNewIntValue(newValue));
  } else if( command == binCmd ) {
    h->SetNumberOfBinsE(binCmd->GetNewIntValue(newValue));
  } else if( command == verbCmd ) {
    h->SetVerbose(verbCmd->GetNewIntValue(newValue));
  } else if (command == beamCmd) {
    h->SetDefaultBeamPositionFlag(beamCmd->GetNewBoolValue(newValue));
  } else if (command == edepCmd) {
    h->SetMaxEnergyDeposit(edepCmd->GetNewDoubleValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

