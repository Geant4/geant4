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
/// \file DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger( DetectorConstruction* Det)
: fDetector(Det)
{
  G4bool broadcast = false;
  fDetDir = new G4UIdirectory("/det/",broadcast);
  fDetDir->SetGuidance("Detector control");

  fAbsMaterCmd = new G4UIcmdWithAString("/det/setAbsMat",this);
  fAbsMaterCmd->SetGuidance("Select Material of the Absorber.");
  fAbsMaterCmd->SetParameterName("AbsMat",false);
  fAbsMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fGapMaterCmd = new G4UIcmdWithAString("/det/setGapMat",this);
  fGapMaterCmd->SetGuidance("Select Material of the Gap.");
  fGapMaterCmd->SetParameterName("GapMat",false);
  fGapMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fAbsThickCmd = new G4UIcmdWithADoubleAndUnit("/det/setAbsThick",this);
  fAbsThickCmd->SetGuidance("Set Thickness of the Absorber");
  fAbsThickCmd->SetParameterName("AbsThick",false);
  fAbsThickCmd->SetRange("AbsThick>=0.");
  fAbsThickCmd->SetUnitCategory("Length");
  fAbsThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fGapThickCmd = new G4UIcmdWithADoubleAndUnit("/det/setGapThick",this);
  fGapThickCmd->SetGuidance("Set Thickness of the Gap");
  fGapThickCmd->SetParameterName("GapThick",false);
  fGapThickCmd->SetRange("GapThick>=0.");
  fGapThickCmd->SetUnitCategory("Length");
  fGapThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fSizeYZCmd = new G4UIcmdWithADoubleAndUnit("/det/setSizeYZ",this);
  fSizeYZCmd->SetGuidance("Set tranverse size of the calorimeter");
  fSizeYZCmd->SetParameterName("SizeYZ",false);
  fSizeYZCmd->SetRange("SizeYZ>0.");
  fSizeYZCmd->SetUnitCategory("Length");
  fSizeYZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fNbLayersCmd = new G4UIcmdWithAnInteger("/det/setNbOfLayers",this);
  fNbLayersCmd->SetGuidance("Set number of layers.");
  fNbLayersCmd->SetParameterName("NbLayers",false);
  fNbLayersCmd->SetRange("NbLayers>0 && NbLayers<500");
  fNbLayersCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fNbLayersCmd;
  delete fAbsMaterCmd;
  delete fGapMaterCmd;
  delete fAbsThickCmd;
  delete fGapThickCmd;
  delete fSizeYZCmd;
  delete fDetDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if( command == fAbsMaterCmd ) {
    fDetector->SetAbsorberMaterial(newValue);
    return;
  }

  if( command == fGapMaterCmd ) {
    fDetector->SetGapMaterial(newValue);
    return;
  }

  if( command == fAbsThickCmd ) {
    fDetector->SetAbsorberThickness(fAbsThickCmd->GetNewDoubleValue(newValue));
    return;
  }

  if( command == fGapThickCmd ) {
    fDetector->SetGapThickness(fGapThickCmd->GetNewDoubleValue(newValue));
    return;
  }

  if( command == fSizeYZCmd ) {
    fDetector->SetCalorSizeYZ(fSizeYZCmd->GetNewDoubleValue(newValue));
    return;
  }

  if( command == fNbLayersCmd ) {
    fDetector->SetNbOfLayers(fNbLayersCmd->GetNewIntValue(newValue));
    return;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
