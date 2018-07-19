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
// Please cite the following paper if you use this software
// Nucl.Instrum.Meth.B260:20-27, 2007

#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:fDetector(Det)
{ 
  fQuadDir = new G4UIdirectory("/quad/");
  fQuadDir->SetGuidance("Quadrupole control.");
  
  fModelCmd = new G4UIcmdWithAnInteger("/quad/setModel",this);
  fModelCmd->SetGuidance("Select magnetic field model.");
  fModelCmd->SetParameterName("model",false);
  fModelCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fProfileCmd = new G4UIcmdWithAnInteger("/displayProfile",this);
  fProfileCmd->SetGuidance("Display beam profile.");
  fProfileCmd->SetParameterName("profile",false);
  fProfileCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fGridCmd = new G4UIcmdWithAnInteger("/setGrid",this);
  fGridCmd->SetGuidance("Put grid and shadow plane.");
  fGridCmd->SetParameterName("grid",false);
  fGridCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fCoefCmd = new G4UIcmdWithAnInteger("/setCoef",this);
  fCoefCmd->SetGuidance("Calculate aberration coefficients.");
  fCoefCmd->SetParameterName("coef",false);
  fCoefCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fG1Cmd = new G4UIcmdWithADouble("/quad/setG1",this);
  fG1Cmd->SetGuidance("Set G1.");
  fG1Cmd->SetParameterName("G1",false);
  fG1Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fG2Cmd = new G4UIcmdWithADouble("/quad/setG2",this);
  fG2Cmd->SetGuidance("Set G2.");
  fG2Cmd->SetParameterName("G2",false);
  fG2Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fG3Cmd = new G4UIcmdWithADouble("/quad/setG3",this);
  fG3Cmd->SetGuidance("Set G3.");
  fG3Cmd->SetParameterName("G3",false);
  fG3Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fG4Cmd = new G4UIcmdWithADouble("/quad/setG4",this);
  fG4Cmd->SetGuidance("Set G4.");
  fG4Cmd->SetParameterName("G4",false);
  fG4Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fG1Cmd;
  delete fG2Cmd;
  delete fG3Cmd;
  delete fG4Cmd;
  delete fQuadDir;  
  delete fModelCmd;
  delete fProfileCmd;
  delete fGridCmd;
  delete fCoefCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == fModelCmd )
   { fDetector->SetModel(fModelCmd->GetNewIntValue(newValue));}

  if( command == fProfileCmd )
   { fDetector->SetProfile(fProfileCmd->GetNewIntValue(newValue));}

  if( command == fGridCmd )
   { fDetector->SetGrid(fGridCmd->GetNewIntValue(newValue));}

  if( command == fCoefCmd )
   { fDetector->SetCoef(fGridCmd->GetNewIntValue(newValue));}

  if( command == fG1Cmd )
   { fDetector->SetG1(fG1Cmd->GetNewDoubleValue(newValue));}
   
  if( command == fG2Cmd )
   { fDetector->SetG2(fG1Cmd->GetNewDoubleValue(newValue));}
   
  if( command == fG3Cmd )
   { fDetector->SetG3(fG1Cmd->GetNewDoubleValue(newValue));}
   
  if( command == fG4Cmd )
   { fDetector->SetG4(fG1Cmd->GetNewDoubleValue(newValue));}
   
}
