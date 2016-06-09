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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "exrdmDetectorMessenger.hh"

#include "exrdmDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

exrdmDetectorMessenger::exrdmDetectorMessenger(exrdmDetectorConstruction* myDet)
:myDetector(myDet)
{ 
  exrdmDir = new G4UIdirectory("/exrdm/");
  exrdmDir->SetGuidance("UI commands specific to this example.");
  
  detDir = new G4UIdirectory("/exrdm/det/");
  detDir->SetGuidance("detector control.");
  
  TargMatCmd = new G4UIcmdWithAString("/exrdm/det/setTargetMate",this);
  TargMatCmd->SetGuidance("Select Material of the Target.");
  TargMatCmd->SetParameterName("choice",false);
  TargMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  TargRadiusCmd = new G4UIcmdWithADoubleAndUnit("/exrdm/det/setTargetRadius",this);
  TargRadiusCmd->SetGuidance("Set the Target Radius.");
  TargRadiusCmd->SetUnitCategory("Length");
  TargRadiusCmd->SetParameterName("choice",false);
  TargRadiusCmd->AvailableForStates(G4State_PreInit);
  
  TargLengthCmd = new G4UIcmdWithADoubleAndUnit("/exrdm/det/setTargetLength",this);
  TargLengthCmd->SetGuidance("Set the Target Length.");
  TargLengthCmd->SetUnitCategory("Length");
  TargLengthCmd->SetParameterName("choice",false);
  TargLengthCmd->AvailableForStates(G4State_PreInit);

  DetectMatCmd = new G4UIcmdWithAString("/exrdm/det/setDetectorMate",this);
  DetectMatCmd->SetGuidance("Select Material of the Detector.");
  DetectMatCmd->SetParameterName("choice",false);
  DetectMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  DetectThicknessCmd = new G4UIcmdWithADoubleAndUnit("/exrdm/det/setDetectorThickness",this);
  DetectThicknessCmd->SetGuidance("Set the Detector Thickness.");
  DetectThicknessCmd->SetUnitCategory("Length");
  DetectThicknessCmd->SetParameterName("choice",false);
  DetectThicknessCmd->AvailableForStates(G4State_PreInit);

  DetectLengthCmd = new G4UIcmdWithADoubleAndUnit("/exrdm/det/setDetectorLength",this);
  DetectLengthCmd->SetGuidance("Set the Detector Length.");
  DetectLengthCmd->SetUnitCategory("Length");
  DetectLengthCmd->SetParameterName("choice",false);
  DetectLengthCmd->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

exrdmDetectorMessenger::~exrdmDetectorMessenger()
{
  delete TargMatCmd;
  delete DetectMatCmd;
  delete TargRadiusCmd;
  delete DetectThicknessCmd;
  delete TargLengthCmd;
  delete DetectLengthCmd;
  delete detDir;
  delete exrdmDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void exrdmDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == TargMatCmd )
    { myDetector->setTargetMaterial(newValue);}
  else if ( command == TargLengthCmd ) 
    { myDetector->setTargetLength(TargLengthCmd->GetNewDoubleValue(newValue));}
  else if ( command == TargRadiusCmd ) 
    { myDetector->setTargetRadius(TargLengthCmd->GetNewDoubleValue(newValue));}
  else if( command == DetectMatCmd )
    { myDetector->setDetectorMaterial(newValue);} 
  else if (command == DetectLengthCmd ) 
    { myDetector->setDetectorLength(DetectLengthCmd->GetNewDoubleValue(newValue));}
  else if (command == DetectThicknessCmd ) 
    { myDetector->setDetectorThickness(DetectThicknessCmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
