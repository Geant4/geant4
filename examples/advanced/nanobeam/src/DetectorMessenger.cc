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
// -------------------------------------------------------------------
// $Id: DetectorMessenger.cc,v 1.2 2008-01-25 20:49:24 sincerti Exp $
// -------------------------------------------------------------------

#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:Detector(Det)
{ 
  quadDir = new G4UIdirectory("/quad/");
  quadDir->SetGuidance("Quadrupole control.");
  
  modelCmd = new G4UIcmdWithAnInteger("/quad/setModel",this);
  modelCmd->SetGuidance("Select magnetic field model.");
  modelCmd->SetParameterName("model",false);
  modelCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  profileCmd = new G4UIcmdWithAnInteger("/displayProfile",this);
  profileCmd->SetGuidance("Display beam profile.");
  profileCmd->SetParameterName("profile",false);
  profileCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  gridCmd = new G4UIcmdWithAnInteger("/setGrid",this);
  gridCmd->SetGuidance("Put grid and shadow plane.");
  gridCmd->SetParameterName("grid",false);
  gridCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  G1Cmd = new G4UIcmdWithADouble("/quad/setG1",this);
  G1Cmd->SetGuidance("Set G1.");
  G1Cmd->SetParameterName("G1",false);
  G1Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  G2Cmd = new G4UIcmdWithADouble("/quad/setG2",this);
  G2Cmd->SetGuidance("Set G2.");
  G2Cmd->SetParameterName("G2",false);
  G2Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  G3Cmd = new G4UIcmdWithADouble("/quad/setG3",this);
  G3Cmd->SetGuidance("Set G3.");
  G3Cmd->SetParameterName("G3",false);
  G3Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  G4Cmd = new G4UIcmdWithADouble("/quad/setG4",this);
  G4Cmd->SetGuidance("Set G4.");
  G4Cmd->SetParameterName("G4",false);
  G4Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  UpdateCmd = new G4UIcmdWithoutParameter("/quad/update",this);
  UpdateCmd->SetGuidance("Update calorimeter geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete G1Cmd;
  delete G2Cmd;
  delete G3Cmd;
  delete G4Cmd;
  delete UpdateCmd;
  delete quadDir;  
  delete modelCmd;
  delete profileCmd;
  delete gridCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == modelCmd )
   { Detector->SetModel(modelCmd->GetNewIntValue(newValue));}

  if( command == profileCmd )
   { Detector->SetProfile(profileCmd->GetNewIntValue(newValue));}

  if( command == gridCmd )
   { Detector->SetGrid(gridCmd->GetNewIntValue(newValue));}

  if( command == G1Cmd )
   { Detector->SetG1(G1Cmd->GetNewDoubleValue(newValue));}
   
  if( command == G2Cmd )
   { Detector->SetG2(G1Cmd->GetNewDoubleValue(newValue));}
   
  if( command == G3Cmd )
   { Detector->SetG3(G1Cmd->GetNewDoubleValue(newValue));}
   
  if( command == G4Cmd )
   { Detector->SetG4(G1Cmd->GetNewDoubleValue(newValue));}
   
  if( command == UpdateCmd )
   { Detector->UpdateGeometry();}
}
