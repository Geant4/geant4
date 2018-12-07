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
/// \file electromagnetic/TestEm8/src/DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
//
/////////////////////////////////////////////////////////////////////////
//
// TestEm8: Gaseous detector
//
// Created: 31.08.2010 V.Ivanchenko
//
// Modified:
//
////////////////////////////////////////////////////////////////////////
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * det)
  : G4UImessenger(),fDetector(det)
{ 
  fDetDir = new G4UIdirectory("/testem/");
  fDetDir->SetGuidance("Detector control.");
      
  fGasMaterCmd = new G4UIcmdWithAString("/testem/setGasMat",this);
  fGasMaterCmd->SetGuidance("Select material of the detector.");
  fGasMaterCmd->SetParameterName("gmat",true);
  fGasMaterCmd->SetDefaultValue("Argon");
  fGasMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fWindowMaterCmd = new G4UIcmdWithAString("/testem/setWindowMat",this);
  fWindowMaterCmd->SetGuidance("Select material of the window.");
  fWindowMaterCmd->SetParameterName("wmat",true);
  fWindowMaterCmd->SetDefaultValue("Mylar");
  fWindowMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fWorldMaterCmd = new G4UIcmdWithAString("/testem/setWorldMat",this);
  fWorldMaterCmd->SetGuidance("Select material of the world.");
  fWorldMaterCmd->SetParameterName("worldmat",true);
  fWorldMaterCmd->SetDefaultValue("empty");
  fWorldMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fGasThickCmd = new G4UIcmdWithADoubleAndUnit("/testem/setGasThick",this);
  fGasThickCmd->SetGuidance("Set thickness of the detector");
  fGasThickCmd->SetParameterName("SizeZ",false,false);
  fGasThickCmd->SetUnitCategory("Length");
  fGasThickCmd->SetDefaultUnit("mm");
  fGasThickCmd->SetRange("SizeZ>0.");
  fGasThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fGasRadCmd = new G4UIcmdWithADoubleAndUnit("/testem/setGasRad",this);
  fGasRadCmd->SetGuidance("Set radius of the detector");
  fGasRadCmd->SetParameterName("SizeR",false,false);
  fGasRadCmd->SetUnitCategory("Length");
  fGasRadCmd->SetDefaultUnit("mm");
  fGasRadCmd->SetRange("SizeR>0.");
  fGasRadCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fWinThickCmd = new G4UIcmdWithADoubleAndUnit("/testem/setWindowThick",this);
  fWinThickCmd->SetGuidance("Set thickness of the window");
  fWinThickCmd->SetParameterName("delta",false,false);
  fWinThickCmd->SetUnitCategory("Length");
  fWinThickCmd->SetDefaultUnit("mm");
  fWinThickCmd->SetRange("delta>0.");
  fWinThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fIonCmd = new G4UIcmdWithADoubleAndUnit("/testem/setPairEnergy",this);
  fIonCmd->SetGuidance("Set energy per electron-ion pair for detector");
  fIonCmd->SetParameterName("en",false,false);
  fIonCmd->SetUnitCategory("Energy");
  fIonCmd->SetDefaultUnit("MeV");
  fIonCmd->SetRange("en>0.");
  fIonCmd->AvailableForStates(G4State_PreInit,G4State_Idle);    

  fStepMaxCmd = new G4UIcmdWithADoubleAndUnit("/testem/stepMax",this);
  fStepMaxCmd->SetGuidance("Set max allowed step length for charged particles");
  fStepMaxCmd->SetParameterName("mxStep",false);
  fStepMaxCmd->SetRange("mxStep>0.");
  fStepMaxCmd->SetUnitCategory("Length");
  fStepMaxCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fGasMaterCmd; 
  delete fGasThickCmd; 
  delete fGasRadCmd;  
  delete fWinThickCmd; 
  delete fWindowMaterCmd;
  delete fWorldMaterCmd;
  delete fIonCmd;
  delete fStepMaxCmd; 
  delete fDetDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == fGasMaterCmd ) { 
    fDetector->SetGasMaterial(newValue);
  } else if( command == fWindowMaterCmd ) { 
    fDetector->SetContainerMaterial(newValue);
  } else if( command == fWorldMaterCmd ) { 
    fDetector->SetWorldMaterial(newValue);
  } else if( command == fGasThickCmd ) { 
    fDetector->SetGasThickness(fGasThickCmd->GetNewDoubleValue(newValue));
  } else if( command == fGasRadCmd ) { 
    fDetector->SetGasRadius(fGasRadCmd->GetNewDoubleValue(newValue));
  } else if( command == fWinThickCmd ) { 
    fDetector->SetContainerThickness(fWinThickCmd->GetNewDoubleValue(newValue));
  } else if( command == fStepMaxCmd ) { 
    fDetector->SetMaxChargedStep(fStepMaxCmd->GetNewDoubleValue(newValue));
  } else if( command == fIonCmd ) { 
    fDetector->SetPairEnergy(fIonCmd->GetNewDoubleValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
