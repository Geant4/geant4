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
#include "Par03DetectorMessenger.hh"
#include "Par03DetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par03DetectorMessenger::Par03DetectorMessenger(
  Par03DetectorConstruction* aDetector)
  : G4UImessenger()
  , fDetector(aDetector)
{
  fExampleDir = new G4UIdirectory("/Par03/");
  fExampleDir->SetGuidance("UI commands specific to this example");

  fDetectorDir = new G4UIdirectory("/Par03/detector/");
  fDetectorDir->SetGuidance("Detector construction UI commands");

  fPrintCmd = new G4UIcmdWithoutParameter("/Par03/detector/print", this);
  fPrintCmd->SetGuidance("Print current settings.");

  fDetectorRadiusCmd =
    new G4UIcmdWithADoubleAndUnit("/Par03/detector/setDetectorRadius", this);
  fDetectorRadiusCmd->SetGuidance(
    "Set tranverse size of the detector (cylinder radius)");
  fDetectorRadiusCmd->SetParameterName("Size", false);
  fDetectorRadiusCmd->SetRange("Size>0.");
  fDetectorRadiusCmd->SetUnitCategory("Length");
  fDetectorRadiusCmd->AvailableForStates(G4State_PreInit);
  fDetectorRadiusCmd->SetToBeBroadcasted(false);

  fDetectorLengthCmd =
    new G4UIcmdWithADoubleAndUnit("/Par03/detector/setDetectorLength", this);
  fDetectorLengthCmd->SetGuidance(
    "Set length of the detector (cylinder length)");
  fDetectorLengthCmd->SetParameterName("Size", false);
  fDetectorLengthCmd->SetRange("Size>0.");
  fDetectorLengthCmd->SetUnitCategory("Length");
  fDetectorLengthCmd->AvailableForStates(G4State_PreInit);
  fDetectorLengthCmd->SetToBeBroadcasted(false);

  fDetectorMaterialCmd =
    new G4UIcmdWithAString("/Par03/detector/setDetectorMaterial", this);
  fDetectorMaterialCmd->SetGuidance("Material of the detector.");
  fDetectorMaterialCmd->SetParameterName("Name", false);
  fDetectorMaterialCmd->AvailableForStates(G4State_PreInit);
  fDetectorMaterialCmd->SetToBeBroadcasted(false);

  fNbLayersCmd =
    new G4UIcmdWithAnInteger("/Par03/detector/setNbOfLayers", this);
  fNbLayersCmd->SetGuidance("Set number of layers.");
  fNbLayersCmd->SetParameterName("NbLayers", false);
  fNbLayersCmd->SetRange("NbLayers>0");
  fNbLayersCmd->AvailableForStates(G4State_PreInit);
  fNbLayersCmd->SetToBeBroadcasted(false);

  fNbRhoCellsCmd =
    new G4UIcmdWithAnInteger("/Par03/detector/setNbOfRhoCells", this);
  fNbRhoCellsCmd->SetGuidance("Set number of cells along radius.");
  fNbRhoCellsCmd->SetParameterName("NbRhoCells", false);
  fNbRhoCellsCmd->SetRange("NbRhoCells>0");
  fNbRhoCellsCmd->AvailableForStates(G4State_PreInit);
  fNbRhoCellsCmd->SetToBeBroadcasted(false);

  fNbPhiCellsCmd =
    new G4UIcmdWithAnInteger("/Par03/detector/setNbOfPhiCells", this);
  fNbPhiCellsCmd->SetGuidance("Set number of cells in azimuthal angle.");
  fNbPhiCellsCmd->SetParameterName("NbPhiCells", false);
  fNbPhiCellsCmd->SetRange("NbPhiCells>0");
  fNbPhiCellsCmd->AvailableForStates(G4State_PreInit);
  fNbPhiCellsCmd->SetToBeBroadcasted(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par03DetectorMessenger::~Par03DetectorMessenger()
{
  delete fPrintCmd;
  delete fDetectorRadiusCmd;
  delete fDetectorLengthCmd;
  delete fDetectorMaterialCmd;
  delete fNbLayersCmd;
  delete fNbRhoCellsCmd;
  delete fNbPhiCellsCmd;
  delete fDetectorDir;
  delete fExampleDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par03DetectorMessenger::SetNewValue(G4UIcommand* aCommand,
                                         G4String aNewValue)
{
  if(aCommand == fPrintCmd)
  {
    fDetector->Print();
  }
  else if(aCommand == fDetectorRadiusCmd)
  {
    fDetector->SetRadius(fDetectorRadiusCmd->GetNewDoubleValue(aNewValue));
  }
  else if(aCommand == fDetectorLengthCmd)
  {
    fDetector->SetLength(fDetectorRadiusCmd->GetNewDoubleValue(aNewValue));
  }
  else if(aCommand == fDetectorMaterialCmd)
  {
    fDetector->SetMaterial(aNewValue);
  }
  else if(aCommand == fNbLayersCmd)
  {
    fDetector->SetNbOfLayers(fNbLayersCmd->GetNewIntValue(aNewValue));
  }
  else if(aCommand == fNbRhoCellsCmd)
  {
    fDetector->SetNbOfRhoCells(fNbRhoCellsCmd->GetNewIntValue(aNewValue));
  }
  else if(aCommand == fNbPhiCellsCmd)
  {
    fDetector->SetNbOfPhiCells(fNbPhiCellsCmd->GetNewIntValue(aNewValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String Par03DetectorMessenger::GetCurrentValue(G4UIcommand* aCommand)
{
  G4String cv;

  if(aCommand == fDetectorRadiusCmd)
  {
    cv = fDetectorRadiusCmd->ConvertToString(fDetector->GetRadius(), "mm");
  }
  else if(aCommand == fDetectorLengthCmd)
  {
    cv = fDetectorLengthCmd->ConvertToString(fDetector->GetLength(), "mm");
  }
  else if(aCommand == fDetectorMaterialCmd)
  {
    cv = fDetector->GetMaterial();
  }
  else if(aCommand == fNbLayersCmd)
  {
    cv = fNbLayersCmd->ConvertToString(fDetector->GetNbOfLayers());
  }
  else if(aCommand == fNbPhiCellsCmd)
  {
    cv = fNbPhiCellsCmd->ConvertToString(fDetector->GetNbOfPhiCells());
  }
  else if(aCommand == fNbRhoCellsCmd)
  {
    cv = fNbRhoCellsCmd->ConvertToString(fDetector->GetNbOfRhoCells());
  }
  return cv;
}