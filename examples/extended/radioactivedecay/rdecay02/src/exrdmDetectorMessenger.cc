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
/// \file radioactivedecay/rdecay02/src/exrdmDetectorMessenger.cc
/// \brief Implementation of the exrdmDetectorMessenger class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "exrdmDetectorMessenger.hh"

#include "exrdmDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

exrdmDetectorMessenger::exrdmDetectorMessenger(exrdmDetectorConstruction* myDet)
:fMyDetector(myDet)
{ 
  fExrdmDir = new G4UIdirectory("/exrdm/");
  fExrdmDir->SetGuidance("UI commands specific to this example.");
  
  fDetDir = new G4UIdirectory("/exrdm/det/");
  fDetDir->SetGuidance("detector control.");
  
  fTargMatCmd = new G4UIcmdWithAString("/exrdm/det/setTargetMate",this);
  fTargMatCmd->SetGuidance("Select Material of the Target.");
  fTargMatCmd->SetParameterName("choice",false);
  fTargMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fTargRadiusCmd = new G4UIcmdWithADoubleAndUnit("/exrdm/det/setTargetRadius", this);
  fTargRadiusCmd->SetGuidance("Set the Target Radius.");
  fTargRadiusCmd->SetUnitCategory("Length");
  fTargRadiusCmd->SetParameterName("choice",false);
  fTargRadiusCmd->AvailableForStates(G4State_PreInit);
  
  fTargLengthCmd = new G4UIcmdWithADoubleAndUnit("/exrdm/det/setTargetLength", this);
  fTargLengthCmd->SetGuidance("Set the Target Length.");
  fTargLengthCmd->SetUnitCategory("Length");
  fTargLengthCmd->SetParameterName("choice",false);
  fTargLengthCmd->AvailableForStates(G4State_PreInit);

  fDetectMatCmd = new G4UIcmdWithAString("/exrdm/det/setDetectorMate",this);
  fDetectMatCmd->SetGuidance("Select Material of the Detector.");
  fDetectMatCmd->SetParameterName("choice",false);
  fDetectMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  fDetectThicknessCmd =
                 new G4UIcmdWithADoubleAndUnit("/exrdm/det/setDetectorThickness",this);
  fDetectThicknessCmd->SetGuidance("Set the Detector Thickness.");
  fDetectThicknessCmd->SetUnitCategory("Length");
  fDetectThicknessCmd->SetParameterName("choice",false);
  fDetectThicknessCmd->AvailableForStates(G4State_PreInit);

  fDetectLengthCmd =
                 new G4UIcmdWithADoubleAndUnit("/exrdm/det/setDetectorLength",this);
  fDetectLengthCmd->SetGuidance("Set the Detector Length.");
  fDetectLengthCmd->SetUnitCategory("Length");
  fDetectLengthCmd->SetParameterName("choice",false);
  fDetectLengthCmd->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

exrdmDetectorMessenger::~exrdmDetectorMessenger()
{
  delete fTargMatCmd;
  delete fDetectMatCmd;
  delete fTargRadiusCmd;
  delete fDetectThicknessCmd;
  delete fTargLengthCmd;
  delete fDetectLengthCmd;
  delete fDetDir;
  delete fExrdmDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void exrdmDetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if( command == fTargMatCmd )fMyDetector->SetTargetMaterial(newValue);
  else if ( command == fTargLengthCmd ) 
    fMyDetector->SetTargetLength(fTargLengthCmd->GetNewDoubleValue(newValue));
  else if ( command == fTargRadiusCmd ) 
    fMyDetector->SetTargetRadius(fTargLengthCmd->GetNewDoubleValue(newValue));
  else if( command == fDetectMatCmd )
    fMyDetector->SetDetectorMaterial(newValue);
  else if (command == fDetectLengthCmd ) 
    fMyDetector->SetDetectorLength(
                     fDetectLengthCmd->GetNewDoubleValue(newValue));
  else if (command == fDetectThicknessCmd ) 
    fMyDetector->SetDetectorThickness(
                              fDetectThicknessCmd->GetNewDoubleValue(newValue));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
